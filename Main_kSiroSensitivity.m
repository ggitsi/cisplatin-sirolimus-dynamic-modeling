function Main_kSiroSensitivity
% STEP 4B: Sensitivity analysis of the uncertain sirolimus coefficient.
% Required file: Fun_Sirolimus.m
%
% This script varies both kSiro and sirolimus concentration while keeping
% cisplatin at zero. It tests whether the modeled direction of effect is
% robust and quantifies how strongly the predicted magnitude depends on
% the assumed kSiro value.

clc;
close all;

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put Main_kSiroSensitivity.m ' ...
           'and Fun_Sirolimus.m in the same MATLAB folder.']);
end

%% Sensitivity settings
C = 0;
kSiroValues = [0.0025 0.005 0.01 0.0255]; % exploratory range, 1/uM
sirolimusDoses = [0 2 5 10 15 20];         % uM
evaluationTimes = [24 48];                  % hours
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

%% Fixed population used for every sensitivity condition
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1+(sdBCL2/meanBCL2)^2));
muLog = log(meanBCL2)-0.5*sigmaLog^2;

rng(randomSeed,'twister');
BCLmitValues = exp(muLog+sigmaLog*randn(nCells,1));

baselineCells = BCLmitValues <= y0(10);
nBaseline = sum(baselineCells);
baselinePercent = 100*nBaseline/nCells;

%% Allocate outputs: kSiro x sirolimus dose x evaluation time
nK = numel(kSiroValues);
nDoses = numel(sirolimusDoses);
nTimes = numel(evaluationTimes);

meanMTOR = zeros(nK,nDoses,nTimes);
meanAutophagy = zeros(nK,nDoses,nTimes);
apoptosis = zeros(nK,nDoses,nTimes);
autophagyFold = zeros(nK,nDoses,nTimes);
controlCorrectedApoptosis = zeros(nK,nDoses,nTimes);

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% Run sensitivity grid
fprintf('============================================================\n');
fprintf('STEP 4B - kSiro SENSITIVITY ANALYSIS\n');
fprintf('Cisplatin C           : %.1f\n',C);
fprintf('Cells                 : %d\n',nCells);
fprintf('Random seed           : %d\n',randomSeed);
fprintf('kSiro values          :');
fprintf(' %.4f',kSiroValues);
fprintf('\n');
fprintf('Sirolimus doses (uM)  :');
fprintf(' %.1f',sirolimusDoses);
fprintf('\n');
fprintf('Total ODE simulations : %d\n',nK*nDoses*nCells);
fprintf('Baseline apoptosis    : %d/%d (%.1f %%)\n', ...
    nBaseline,nCells,baselinePercent);
fprintf('============================================================\n');

for iK = 1:nK
    kSiro = kSiroValues(iK);

    for iDose = 1:nDoses
        Sirolimus = sirolimusDoses(iDose);
        deathTime = nan(nCells,1);
        mTORAtEval = zeros(nCells,nTimes);
        autophagyAtEval = zeros(nCells,nTimes);

        fprintf('Running kSiro = %.4f, sirolimus = %.1f uM ...\n', ...
            kSiro,Sirolimus);

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
                timeIndex = find(abs(tCell-evaluationTimes(iTime)) < 1e-9, ...
                    1,'first');
                mTORAtEval(iCell,iTime) = yCell(timeIndex,7);
                autophagyAtEval(iCell,iTime) = yCell(timeIndex,2);
            end
        end

        for iTime = 1:nTimes
            meanMTOR(iK,iDose,iTime) = mean(mTORAtEval(:,iTime));
            meanAutophagy(iK,iDose,iTime) = ...
                mean(autophagyAtEval(:,iTime));
            apoptosis(iK,iDose,iTime) = ...
                100*mean(deathTime <= evaluationTimes(iTime));
        end
    end
end

%% Normalize within each kSiro scenario using sirolimus = 0
for iK = 1:nK
    for iTime = 1:nTimes
        controlAutophagy = meanAutophagy(iK,1,iTime);
        autophagyFold(iK,:,iTime) = ...
            meanAutophagy(iK,:,iTime)/controlAutophagy;
        controlCorrectedApoptosis(iK,:,iTime) = ...
            apoptosis(iK,:,iTime)-apoptosis(iK,1,iTime);
    end
end

%% Long-format results table
[kGrid,doseGrid] = ndgrid(kSiroValues,sirolimusDoses);

results = table( ...
    kGrid(:),doseGrid(:), ...
    reshape(meanMTOR(:,:,1),[],1),reshape(meanMTOR(:,:,2),[],1), ...
    reshape(meanAutophagy(:,:,1),[],1), ...
    reshape(meanAutophagy(:,:,2),[],1), ...
    reshape(autophagyFold(:,:,1),[],1), ...
    reshape(autophagyFold(:,:,2),[],1), ...
    reshape(apoptosis(:,:,1),[],1),reshape(apoptosis(:,:,2),[],1), ...
    reshape(controlCorrectedApoptosis(:,:,1),[],1), ...
    reshape(controlCorrectedApoptosis(:,:,2),[],1), ...
    'VariableNames',{ ...
    'kSiro_per_uM','Sirolimus_uM', ...
    'MeanMTOR_24h','MeanMTOR_48h', ...
    'MeanAutophagy_24h','MeanAutophagy_48h', ...
    'AutophagyFold_24h','AutophagyFold_48h', ...
    'Apoptosis_24h_pct','Apoptosis_48h_pct', ...
    'ControlCorrected_24h_pctPoint', ...
    'ControlCorrected_48h_pctPoint'});

fprintf('\nFULL SENSITIVITY RESULTS\n');
disp(results);

%% Summary at the highest sirolimus concentration
summary20uM = table( ...
    kSiroValues.', ...
    meanMTOR(:,end,2), ...
    meanAutophagy(:,end,2), ...
    autophagyFold(:,end,2), ...
    apoptosis(:,end,2), ...
    'VariableNames',{ ...
    'kSiro_per_uM','MeanMTOR_48h','MeanAutophagy_48h', ...
    'AutophagyFold_48h','Apoptosis_48h_pct'});

fprintf('\nSUMMARY AT SIROLIMUS = 20 uM, TIME = 48 h\n');
disp(summary20uM);

directionPass = true(nK,1);
for iK = 1:nK
    directionPass(iK) = ...
        meanMTOR(iK,end,2) < meanMTOR(iK,1,2) && ...
        meanAutophagy(iK,end,2) > meanAutophagy(iK,1,2);
end

if all(directionPass)
    fprintf('Direction robust across all kSiro values: YES\n');
else
    fprintf('Direction robust across all kSiro values: NO\n');
end

outputFile = 'Step4B_kSiro_Sensitivity_Results.csv';
writetable(results,outputFile);
fprintf('Results saved to: %s\n',outputFile);
fprintf('============================================================\n');

%% Sensitivity figures at 48 hours
fig = figure('Color','w','Name','Step 4B - kSiro Sensitivity');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

grayLevel = linspace(0,0.65,nK);
lineStyles = {'-','--',':','-.'};
markers = {'o','s','d','^'};

nexttile;
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,meanMTOR(iK,:,2), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
xlabel('Sirolimus (uM)');
ylabel('Mean mTOR at 48 h');
title('mTOR Sensitivity');
legend('Location','best');
grid on;

nexttile;
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,meanAutophagy(iK,:,2), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
xlabel('Sirolimus (uM)');
ylabel('Mean ATPHG at 48 h');
title('Autophagy Magnitude Sensitivity');
legend('Location','best');
grid on;

nexttile;
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,autophagyFold(iK,:,2), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
xlabel('Sirolimus (uM)');
ylabel('Autophagy fold-change at 48 h');
title('Autophagy Relative to 0 uM');
legend('Location','best');
grid on;

nexttile;
hold on;
for iK = 1:nK
    curveColor = grayLevel(iK)*[1 1 1];
    plot(sirolimusDoses,apoptosis(iK,:,2), ...
        'Color',curveColor,'LineStyle',lineStyles{iK}, ...
        'Marker',markers{iK},'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',sprintf('k = %.4f',kSiroValues(iK)));
end
xlabel('Sirolimus (uM)');
ylabel('Apoptotic cells at 48 h (%)');
title('Sirolimus-Only Apoptosis');
legend('Location','best');
ylim([0 100]);
grid on;

% Explicit light theme.
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
