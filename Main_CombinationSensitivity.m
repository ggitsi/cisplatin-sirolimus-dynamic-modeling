function Main_CombinationSensitivity
% STEP 5: Cisplatin-sirolimus combination sensitivity analysis.
% Required file: Fun_Sirolimus.m
%
% Positive Protection values mean that adding sirolimus reduces apoptosis
% relative to cisplatin monotherapy at the same cisplatin dose and time.

clc;
close all;

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put Main_CombinationSensitivity.m ' ...
           'and Fun_Sirolimus.m in the same MATLAB folder.']);
end

%% Combination grid
cisplatinDoses = [0 5 6 10 20];       % model input; C=20 matches 20 uM benchmark
sirolimusDoses = [0 2 5 10 20];       % uM
kSiroValues = [0.0025 0.005 0.01 0.0255]; % exploratory range, 1/uM
kSiroNominal = 0.005;                  % conservative nominal scenario
evaluationTimes = [24 48];             % hours
simulationTime = max(evaluationTimes);
nCells = 100;
randomSeed = 42;
timeStep = 0.2;
tEval = 0:timeStep:simulationTime;

%% Fixed parameters from the Tavassoly model
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

%% Fixed population shared by every condition
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1+(sdBCL2/meanBCL2)^2));
muLog = log(meanBCL2)-0.5*sigmaLog^2;

rng(randomSeed,'twister');
BCLmitValues = exp(muLog+sigmaLog*randn(nCells,1));

baselineCells = BCLmitValues <= y0(10);
nBaseline = sum(baselineCells);
baselinePercent = 100*nBaseline/nCells;

%% Allocate outputs: kSiro x cisplatin x sirolimus x time
nK = numel(kSiroValues);
nC = numel(cisplatinDoses);
nS = numel(sirolimusDoses);
nTimes = numel(evaluationTimes);

apoptosis = zeros(nK,nC,nS,nTimes);
meanAutophagy = zeros(nK,nC,nS,nTimes);
meanMTOR = zeros(nK,nC,nS,nTimes);
protection = zeros(nK,nC,nS,nTimes);
netApoptosis = zeros(nK,nC,nS,nTimes);
preventedFraction = nan(nK,nC,nS,nTimes);

evaluationIndices = zeros(1,nTimes);
for iTime = 1:nTimes
    evaluationIndices(iTime) = find( ...
        abs(tEval-evaluationTimes(iTime)) < 1e-9,1,'first');
end

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Sirolimus=0 does not depend on kSiro, so it is simulated once and copied.
nUniqueConditions = nC*nS+(nK-1)*nC*(nS-1);
totalODESimulations = nUniqueConditions*nCells;
conditionCounter = 0;

fprintf('============================================================\n');
fprintf('STEP 5 - CISPLATIN-SIROLIMUS COMBINATION SENSITIVITY\n');
fprintf('Cells                  : %d\n',nCells);
fprintf('Random seed            : %d\n',randomSeed);
fprintf('Cisplatin doses        :');
fprintf(' %.1f',cisplatinDoses);
fprintf('\n');
fprintf('Sirolimus doses (uM)   :');
fprintf(' %.1f',sirolimusDoses);
fprintf('\n');
fprintf('kSiro values           :');
fprintf(' %.4f',kSiroValues);
fprintf('\n');
fprintf('Unique conditions      : %d\n',nUniqueConditions);
fprintf('Total ODE simulations  : %d\n',totalODESimulations);
fprintf('Baseline apoptosis     : %d/%d (%.1f %%)\n', ...
    nBaseline,nCells,baselinePercent);
fprintf('============================================================\n');

%% Run the combination grid
for iK = 1:nK
    kSiro = kSiroValues(iK);

    for iC = 1:nC
        C = cisplatinDoses(iC);

        for iS = 1:nS
            Sirolimus = sirolimusDoses(iS);

            % Reuse monotherapy because kSiro*0 = 0 for every kSiro.
            if iK > 1 && iS == 1
                apoptosis(iK,iC,iS,:) = apoptosis(1,iC,iS,:);
                meanAutophagy(iK,iC,iS,:) = meanAutophagy(1,iC,iS,:);
                meanMTOR(iK,iC,iS,:) = meanMTOR(1,iC,iS,:);
                continue;
            end

            conditionCounter = conditionCounter+1;
            fprintf(['Condition %d/%d: kSiro=%.4f, ' ...
                     'cisplatin=%.1f, sirolimus=%.1f uM\n'], ...
                conditionCounter,nUniqueConditions,kSiro,C,Sirolimus);

            deathTime = nan(nCells,1);
            autophagyAtEval = zeros(nCells,nTimes);
            mTORAtEval = zeros(nCells,nTimes);

            for iCell = 1:nCells
                BCLmit = BCLmitValues(iCell);

                [tCell,yCell] = ode15s(@Fun_Sirolimus,tEval,y0,options, ...
                    ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl, ...
                    sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal, ...
                    kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0, ...
                    wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m, ...
                    wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP, ...
                    Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp, ...
                    sigmaA,wbh_c,wbh_a,Sirolimus,kSiro);

                firstCrossing = find(yCell(:,10) >= BCLmit,1,'first');
                if ~isempty(firstCrossing)
                    deathTime(iCell) = tCell(firstCrossing);
                end

                for iTime = 1:nTimes
                    idx = evaluationIndices(iTime);
                    autophagyAtEval(iCell,iTime) = yCell(idx,2);
                    mTORAtEval(iCell,iTime) = yCell(idx,7);
                end
            end

            for iTime = 1:nTimes
                apoptosis(iK,iC,iS,iTime) = ...
                    100*mean(deathTime <= evaluationTimes(iTime));
                meanAutophagy(iK,iC,iS,iTime) = ...
                    mean(autophagyAtEval(:,iTime));
                meanMTOR(iK,iC,iS,iTime) = ...
                    mean(mTORAtEval(:,iTime));
            end
        end
    end
end

%% Calculate treatment effects
for iK = 1:nK
    for iTime = 1:nTimes
        baselineControl = apoptosis(iK,1,1,iTime);

        for iC = 1:nC
            monoApoptosis = apoptosis(iK,iC,1,iTime);
            inducedByCisplatin = monoApoptosis-baselineControl;

            for iS = 1:nS
                comboApoptosis = apoptosis(iK,iC,iS,iTime);

                % Positive means protection; negative means more apoptosis.
                protection(iK,iC,iS,iTime) = ...
                    monoApoptosis-comboApoptosis;

                netApoptosis(iK,iC,iS,iTime) = ...
                    comboApoptosis-baselineControl;

                if inducedByCisplatin > 0
                    preventedFraction(iK,iC,iS,iTime) = ...
                        100*protection(iK,iC,iS,iTime)/inducedByCisplatin;
                end
            end
        end
    end
end

%% Long-format table
[kGrid,cGrid,sGrid] = ndgrid( ...
    kSiroValues,cisplatinDoses,sirolimusDoses);

results = table( ...
    kGrid(:),cGrid(:),sGrid(:), ...
    reshape(apoptosis(:,:,:,1),[],1), ...
    reshape(apoptosis(:,:,:,2),[],1), ...
    reshape(netApoptosis(:,:,:,1),[],1), ...
    reshape(netApoptosis(:,:,:,2),[],1), ...
    reshape(protection(:,:,:,1),[],1), ...
    reshape(protection(:,:,:,2),[],1), ...
    reshape(preventedFraction(:,:,:,1),[],1), ...
    reshape(preventedFraction(:,:,:,2),[],1), ...
    reshape(meanAutophagy(:,:,:,1),[],1), ...
    reshape(meanAutophagy(:,:,:,2),[],1), ...
    reshape(meanMTOR(:,:,:,1),[],1), ...
    reshape(meanMTOR(:,:,:,2),[],1), ...
    'VariableNames',{ ...
    'kSiro_per_uM','Cisplatin_C','Sirolimus_uM', ...
    'Apoptosis_24h_pct','Apoptosis_48h_pct', ...
    'NetApoptosis_24h_pctPoint','NetApoptosis_48h_pctPoint', ...
    'Protection_24h_pctPoint','Protection_48h_pctPoint', ...
    'PreventedCisplatinEffect_24h_pct', ...
    'PreventedCisplatinEffect_48h_pct', ...
    'MeanAutophagy_24h','MeanAutophagy_48h', ...
    'MeanMTOR_24h','MeanMTOR_48h'});

outputFile = 'Step5_CombinationSensitivity_Results.csv';
writetable(results,outputFile);

%% Robustness table across all kSiro values at 48 hours
nRobustRows = (nC-1)*(nS-1);
robustC = zeros(nRobustRows,1);
robustS = zeros(nRobustRows,1);
minProtection = zeros(nRobustRows,1);
maxProtection = zeros(nRobustRows,1);
robustlyProtective = false(nRobustRows,1);

iRow = 0;
for iC = 2:nC
    for iS = 2:nS
        iRow = iRow+1;
        values = reshape(protection(:,iC,iS,2),[],1);
        robustC(iRow) = cisplatinDoses(iC);
        robustS(iRow) = sirolimusDoses(iS);
        minProtection(iRow) = min(values);
        maxProtection(iRow) = max(values);
        robustlyProtective(iRow) = all(values > 0);
    end
end

robustness = table( ...
    robustC,robustS,minProtection,maxProtection,robustlyProtective, ...
    'VariableNames',{ ...
    'Cisplatin_C','Sirolimus_uM', ...
    'MinProtection_48h_pctPoint','MaxProtection_48h_pctPoint', ...
    'ProtectiveForAll_kSiro'});

robustnessFile = 'Step5_RobustnessSummary_48h.csv';
writetable(robustness,robustnessFile);

fprintf('\nROBUSTNESS SUMMARY AT 48 HOURS\n');
disp(robustness);

nRobust = sum(robustlyProtective);
fprintf('Robustly protective combinations: %d/%d\n', ...
    nRobust,nRobustRows);

%% Focused summary at high cisplatin and 20 uM sirolimus
focusC = [10 20];
nFocus = nK*numel(focusC);
focusK = zeros(nFocus,1);
focusCOut = zeros(nFocus,1);
focusMono = zeros(nFocus,1);
focusCombo = zeros(nFocus,1);
focusProtection = zeros(nFocus,1);
focusPrevented = zeros(nFocus,1);

iRow = 0;
iS20 = find(sirolimusDoses == 20,1);
for iK = 1:nK
    for iFocus = 1:numel(focusC)
        iRow = iRow+1;
        iC = find(cisplatinDoses == focusC(iFocus),1);
        focusK(iRow) = kSiroValues(iK);
        focusCOut(iRow) = focusC(iFocus);
        focusMono(iRow) = apoptosis(iK,iC,1,2);
        focusCombo(iRow) = apoptosis(iK,iC,iS20,2);
        focusProtection(iRow) = protection(iK,iC,iS20,2);
        focusPrevented(iRow) = preventedFraction(iK,iC,iS20,2);
    end
end

focusedSummary = table( ...
    focusK,focusCOut,focusMono,focusCombo,focusProtection,focusPrevented, ...
    'VariableNames',{ ...
    'kSiro_per_uM','Cisplatin_C','MonotherapyApoptosis_48h_pct', ...
    'CombinationApoptosis_48h_pct','Protection_48h_pctPoint', ...
    'PreventedCisplatinEffect_48h_pct'});

fprintf('\nFOCUSED SUMMARY: SIROLIMUS = 20 uM, TIME = 48 h\n');
disp(focusedSummary);

save('Step5_CombinationSensitivity_Workspace.mat', ...
    'results','robustness','focusedSummary', ...
    'apoptosis','netApoptosis','protection','preventedFraction', ...
    'meanAutophagy','meanMTOR','BCLmitValues', ...
    'cisplatinDoses','sirolimusDoses','kSiroValues','evaluationTimes');

fprintf('Full results saved to : %s\n',outputFile);
fprintf('Robustness saved to   : %s\n',robustnessFile);
fprintf('Workspace saved to    : Step5_CombinationSensitivity_Workspace.mat\n');
fprintf('============================================================\n');

%% Figures
iNominal = find(abs(kSiroValues-kSiroNominal) < 1e-12,1);

fig = figure('Color','w','Name','Step 5 - Combination Sensitivity');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
matrix24 = squeeze(apoptosis(iNominal,:,:,1));
imagesc(sirolimusDoses,cisplatinDoses,matrix24);
set(gca,'YDir','normal');
xlabel('Sirolimus (uM)');
ylabel('Cisplatin input, C');
title(sprintf('Apoptosis at 24 h, k = %.4f',kSiroNominal));
colorbar;
caxis([0 100]);

nexttile;
matrix48 = squeeze(apoptosis(iNominal,:,:,2));
imagesc(sirolimusDoses,cisplatinDoses,matrix48);
set(gca,'YDir','normal');
xlabel('Sirolimus (uM)');
ylabel('Cisplatin input, C');
title(sprintf('Apoptosis at 48 h, k = %.4f',kSiroNominal));
colorbar;
caxis([0 100]);

grayLevel = linspace(0,0.65,nK);
lineStyles = {'-','--',':','-.'};
markers = {'o','s','d','^'};

nexttile;
iC10 = find(cisplatinDoses == 10,1);
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,squeeze(protection(iK,iC10,:,2)), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
yline(0,':k');
xlabel('Sirolimus (uM)');
ylabel('Protection at 48 h (percentage points)');
title('Protection at Cisplatin C = 10');
legend('Location','best');
grid on;

nexttile;
iC20 = find(cisplatinDoses == 20,1);
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,squeeze(protection(iK,iC20,:,2)), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
yline(0,':k');
xlabel('Sirolimus (uM)');
ylabel('Protection at 48 h (percentage points)');
title('Protection at Cisplatin C = 20');
legend('Location','best');
grid on;

colormap(fig,flipud(gray(256)));

allAxes = findall(fig,'Type','axes');
set(allAxes, ...
    'Color','w','XColor','k','YColor','k', ...
    'GridColor',[0.70 0.70 0.70], ...
    'MinorGridColor',[0.85 0.85 0.85], ...
    'Box','on');
set(findall(fig,'Type','text'),'Color','k');

allLegends = findall(fig,'Type','legend');
set(allLegends,'Color','w','TextColor','k','EdgeColor','k');

end
