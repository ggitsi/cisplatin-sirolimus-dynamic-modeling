function Main_kSiroOptimizationRobustness
% STEP 8 - ROBUSTNESS OF OPTIMIZATION ACROSS kSiro
% Tests whether Step 7 optimization remains robust across uncertain kSiro.
% Primary endpoint: highest cisplatin model input C with apoptosis <=20% at 48 h.
% Required: Fun_Sirolimus.m

clc;
close all;

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put this file and ' ...
           'Fun_Sirolimus.m in the same MATLAB folder.']);
end

%% Output folder
outputFolder = 'D:\Kuliah\Transfer Kredit\Final Project Dynamic Modelling';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% Settings
cisplatinGrid = 0:0.5:20;
sirolimusDoses = [0 2 5 10 20];
kSiroValues = [0.0025 0.005 0.01 0.0255];
simulationTime = 48;
targetApoptosis = 20;
nCells = 100;
randomSeed = 42;

nK = numel(kSiroValues);
nC = numel(cisplatinGrid);
nS = numel(sirolimusDoses);

%% Fixed parameters from Tavassoly model
BCL2T=3; IP3RT=1; CaT=2; CASP=0;
wjnk_s=1; wdapk_s=1; wbcl_j=1; wmtor_ca=-1; watg_m=-1;
wbh_a=1; wbec_d=1; Wcp_Ca=1; W_ATG5C=-1;
ka=1.77; kda=0.0948; ksb=2.08; kc=0.510; krb=1.41; kra=3.83;
Gj=1.73; Gd=4.05; Gb=5.21; Gt=1.61; Gg=1.04; Gh=0.0101;
Gl=3.43; Ga=0.5240; Gcp=1.95;
sigmaj=2.42; sigmad=1.01; sigmab=4.57; sigmah=2.89;
sigmag=20.8; sigmat=3.51; sigmal=7.99; sigmacp=32.3;
sigmaA=4.83; kcasp=2.01; kin=9.64; kout=6.31;
wjnk0=-1.22; wdapk0=-2.98; wbcl0=-0.614; wmtor0=0.202;
watg0=0.144; wbh0=-1.26; wbec0=-0.647;
W_ATG50=0.215; Wcp_0=-0.514;
wbh_s=0.0801; wbh_c=0.0003;

y0 = [3 0.285 0.831 0.281 0.103 0.537 0.335 ...
      0.0184 0.397 0.069 0.0382 0.0221 0.717];

%% Same 100-cell population
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1 + (sdBCL2/meanBCL2)^2));
muLog = log(meanBCL2) - 0.5*sigmaLog^2;
rng(randomSeed,'twister');
BCLmitValues = exp(muLog + sigmaLog*randn(nCells,1));
baselineApoptosis = 100*mean(BCLmitValues <= y0(10));

%% Allocate
apoptosis = nan(nK,nC,nS);
optionsBase = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Sirolimus=0 is independent of kSiro, so simulate once and reuse.
nUniqueConditions = nC*nS + (nK-1)*nC*(nS-1);
maxODERuns = nUniqueConditions*nCells;
conditionCounter = 0;

fprintf('============================================================\n');
fprintf('STEP 8 - ROBUSTNESS OF OPTIMIZATION ACROSS kSiro\n');
fprintf('Cells                    : %d\n',nCells);
fprintf('Random seed              : %d\n',randomSeed);
fprintf('Simulation time          : %.1f h\n',simulationTime);
fprintf('Target apoptosis         : %.1f %%\n',targetApoptosis);
fprintf('Cisplatin C range        : %.1f to %.1f\n',min(cisplatinGrid),max(cisplatinGrid));
fprintf('Cisplatin grid step      : %.2f\n',cisplatinGrid(2)-cisplatinGrid(1));
fprintf('Sirolimus doses (uM)     :'); fprintf(' %.1f',sirolimusDoses); fprintf('\n');
fprintf('kSiro values (per uM)    :'); fprintf(' %.4f',kSiroValues); fprintf('\n');
fprintf('Baseline apoptosis       : %.1f %%\n',baselineApoptosis);
fprintf('Unique conditions        : %d\n',nUniqueConditions);
fprintf('Maximum ODE simulations  : %d\n',maxODERuns);
fprintf('============================================================\n');

%% Run grid
for iK = 1:nK
    kSiro = kSiroValues(iK);

    for iC = 1:nC
        C = cisplatinGrid(iC);

        for iS = 1:nS
            Sirolimus = sirolimusDoses(iS);

            if iK > 1 && iS == 1
                apoptosis(iK,iC,iS) = apoptosis(1,iC,iS);
                continue;
            end

            conditionCounter = conditionCounter + 1;
            fprintf(['Condition %d/%d: kSiro=%.4f, C=%.2f, ' ...
                     'Sirolimus=%.1f uM\n'], ...
                     conditionCounter,nUniqueConditions,kSiro,C,Sirolimus);

            dead = false(nCells,1);

            for iCell = 1:nCells
                BCLmit = BCLmitValues(iCell);

                if y0(10) >= BCLmit
                    dead(iCell) = true;
                    continue;
                end

                eventFunction = @(t,y) apoptosisEvent(t,y,BCLmit);
                optionsCell = odeset(optionsBase,'Events',eventFunction);

                odeFunction = @(t,y) Fun_Sirolimus( ...
                    t,y,ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl, ...
                    sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal, ...
                    kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0, ...
                    wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m, ...
                    wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP, ...
                    Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp, ...
                    sigmaA,wbh_c,wbh_a,Sirolimus,kSiro);

                [~,~,te] = ode15s(odeFunction,[0 simulationTime],y0,optionsCell);

                if ~isempty(te)
                    dead(iCell) = true;
                end
            end

            apoptosis(iK,iC,iS) = 100*mean(dead);
        end
    end
end

%% Cmax for each kSiro x sirolimus
Cmax = nan(nK,nS);
reachedUpperBound = false(nK,nS);

for iK = 1:nK
    for iS = 1:nS
        currentCurve = squeeze(apoptosis(iK,:,iS));
        idx = find(currentCurve <= targetApoptosis,1,'last');
        if ~isempty(idx)
            Cmax(iK,iS) = cisplatinGrid(idx);
            reachedUpperBound(iK,iS) = (idx == nC);
        end
    end
end

baselineCmax = Cmax(1,1);
increaseVsNoSiro = Cmax - baselineCmax;

%% Long table
[kGrid,sGrid] = ndgrid(kSiroValues,sirolimusDoses);
longTable = table( ...
    kGrid(:),sGrid(:),repmat(targetApoptosis,numel(kGrid),1), ...
    Cmax(:),increaseVsNoSiro(:),reachedUpperBound(:), ...
    'VariableNames',{ ...
    'kSiro_per_uM','Sirolimus_uM','ApoptosisThreshold_pct', ...
    'MaxC_withinScannedGrid','IncreaseInToleratedC_vs_NoSirolimus', ...
    'ReachedScanUpperBound'});

longFile = fullfile(outputFolder,'Step8_kSiro_Cmax_LongTable.csv');
writetable(longTable,longFile);

%% Robust summary across kSiro
robustMinC = min(Cmax,[],1,'omitnan');
robustMaxC = max(Cmax,[],1,'omitnan');
meanC = mean(Cmax,1,'omitnan');
sdC = std(Cmax,0,1,'omitnan');
rangeWidth = robustMaxC - robustMinC;
conservativeGain = robustMinC - baselineCmax;

improvesAllK = false(1,nS);
for iS = 1:nS
    if iS == 1
        improvesAllK(iS) = true;
    else
        improvesAllK(iS) = all(Cmax(:,iS) > baselineCmax);
    end
end

summaryTable = table( ...
    sirolimusDoses.',robustMinC.',robustMaxC.',meanC.',sdC.', ...
    rangeWidth.',conservativeGain.',improvesAllK.', ...
    'VariableNames',{ ...
    'Sirolimus_uM','RobustMinC_across_kSiro','RobustMaxC_across_kSiro', ...
    'MeanCmax_across_kSiro','SD_Cmax_across_kSiro','CmaxRangeWidth', ...
    'ConservativeGain_vs_NoSirolimus','ImprovesCmax_for_All_kSiro'});

summaryFile = fullfile(outputFolder,'Step8_kSiro_RobustnessSummary.csv');
writetable(summaryTable,summaryFile);

%% Full apoptosis grid table
[kFull,cFull,sFull] = ndgrid(kSiroValues,cisplatinGrid,sirolimusDoses);
apoptosisTable = table( ...
    kFull(:),cFull(:),sFull(:),reshape(apoptosis,[],1), ...
    'VariableNames',{ ...
    'kSiro_per_uM','Cisplatin_C','Sirolimus_uM','Apoptosis_48h_pct'});

apoptosisFile = fullfile(outputFolder,'Step8_ApoptosisGrid_48h.csv');
writetable(apoptosisTable,apoptosisFile);

%% Print summary
fprintf('\n============================================================\n');
fprintf('STEP 8 - Cmax FOR EACH kSiro\n');
fprintf('Target: apoptosis <= %.1f %% at 48 h\n',targetApoptosis);
fprintf('============================================================\n');
disp(longTable);

fprintf('\n============================================================\n');
fprintf('ROBUST SUMMARY ACROSS kSiro\n');
fprintf('Baseline Cmax without sirolimus = %.2f\n',baselineCmax);
fprintf('============================================================\n');
disp(summaryTable);

fprintf('\nInterpretation:\n');
fprintf(['RobustMinC_across_kSiro = conservative Cmax achieved for ' ...
         'every tested kSiro value.\n']);
fprintf(['ConservativeGain_vs_NoSirolimus = minimum improvement over ' ...
         'cisplatin alone across kSiro uncertainty.\n']);
fprintf(['C is a model input and must not automatically be interpreted ' ...
         'as cisplatin concentration in uM.\n']);

%% Save workspace
workspaceFile = fullfile(outputFolder,'Step8_kSiroOptimization_Workspace.mat');
save(workspaceFile, ...
    'cisplatinGrid','sirolimusDoses','kSiroValues','simulationTime', ...
    'targetApoptosis','BCLmitValues','baselineApoptosis','apoptosis', ...
    'Cmax','increaseVsNoSiro','baselineCmax','longTable','summaryTable');

%% Figure
fig = figure('Color','w','Name','Step 8 - kSiro Optimization Robustness');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% Panel 1
nexttile;
hold on;
for iK = 1:nK
    plot(sirolimusDoses,Cmax(iK,:),'-o','LineWidth',1.8,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
yline(baselineCmax,'--',sprintf('No-sirolimus Cmax = %.1f',baselineCmax), ...
    'HandleVisibility','off');
xlabel('Sirolimus (uM)');
ylabel('Maximum C with apoptosis <= 20%');
title('Cmax Sensitivity to kSiro');
grid on;
legend('Location','best');

% Panel 2
nexttile;
imagesc(sirolimusDoses,1:nK,Cmax);
set(gca,'YDir','normal');
yticks(1:nK);
yticklabels(compose('%.4f',kSiroValues));
xlabel('Sirolimus (uM)');
ylabel('kSiro (per uM)');
title('Maximum C at 20% Apoptosis Threshold');
cb = colorbar;
cb.Label.String = 'Maximum C';

% Panel 3
nexttile;
hold on;
iS10 = find(sirolimusDoses == 10,1);
for iK = 1:nK
    plot(cisplatinGrid,squeeze(apoptosis(iK,:,iS10)),'-o', ...
        'LineWidth',1.6,'MarkerSize',3.5, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
yline(targetApoptosis,'--',sprintf('Target = %.0f%%',targetApoptosis), ...
    'HandleVisibility','off');
xlabel('Cisplatin model input, C');
ylabel('Apoptotic cells at 48 h (%)');
title('10 uM Sirolimus Across kSiro');
ylim([0 100]);
grid on;
legend('Location','best');

% Panel 4
nexttile;
yyaxis left
plot(sirolimusDoses,robustMinC,'-o','LineWidth',1.8,'MarkerSize',6);
ylabel('Robust minimum Cmax');
yyaxis right
plot(sirolimusDoses,conservativeGain,'-s','LineWidth',1.8,'MarkerSize',6);
ylabel('Conservative gain vs no sirolimus');
xlabel('Sirolimus (uM)');
title('Worst-Case Robustness Across kSiro');
grid on;

allAxes = findall(fig,'Type','axes');
set(allAxes,'Color','w','XColor','k','YColor','k', ...
    'GridColor',[0.75 0.75 0.75],'Box','on');
set(findall(fig,'Type','text'),'Color','k');

figureFile = fullfile(outputFolder,'Step8_kSiroOptimization_Figure.png');
try
    exportgraphics(fig,figureFile,'Resolution',300);
catch
    fprintf('Figure could not be exported automatically.\n');
end

fprintf('\n============================================================\n');
fprintf('FILES SAVED\n');
fprintf('============================================================\n');
fprintf('Cmax long table   : %s\n',longFile);
fprintf('Robust summary    : %s\n',summaryFile);
fprintf('Apoptosis grid    : %s\n',apoptosisFile);
fprintf('Workspace         : %s\n',workspaceFile);
fprintf('Figure            : %s\n',figureFile);
fprintf('============================================================\n');

end

function [value,isterminal,direction] = apoptosisEvent(~,y,BCLmit)
value = y(10) - BCLmit;
isterminal = 1;
direction = 1;
end
