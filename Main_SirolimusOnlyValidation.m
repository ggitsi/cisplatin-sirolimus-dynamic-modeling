function Main_SirolimusOnlyValidation
% STEP 4: Validate the modeled effect of sirolimus without cisplatin.
% Put this file in the same folder as Fun_Sirolimus.m, then run:
%   Main_SirolimusOnlyValidation
%
% Expected qualitative behavior:
%   sirolimus increases -> mTOR decreases -> autophagy increases.

clc;
close all;

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put both Step 4 files ' ...
           'in the same MATLAB folder.']);
end

%% Simulation settings
C = 0;                              % no cisplatin in Step 4
sirolimusDoses = [0 2 5 10 15 20]; % uM, required exploration range
kSiro = 0.0255;                     % 1/uM, exploratory assumption
evaluationTimes = [24 48];          % hours
simulationTime = max(evaluationTimes);
nCells = 100;
randomSeed = 42;
timeStep = 0.1;
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

%% Generate the same virtual population used in Steps 2 and 3
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1+(sdBCL2/meanBCL2)^2));
muLog = log(meanBCL2)-0.5*sigmaLog^2;

rng(randomSeed,'twister');
BCLmitValues = exp(muLog+sigmaLog*randn(nCells,1));

baselineCells = BCLmitValues <= y0(10);
nBaseline = sum(baselineCells);
baselinePercent = 100*nBaseline/nCells;

%% Allocate outputs
nDoses = numel(sirolimusDoses);
nTimes = numel(evaluationTimes);
meanMTORAtTime = zeros(nDoses,nTimes);
meanAutophagyAtTime = zeros(nDoses,nTimes);
apoptosisAtTime = zeros(nDoses,nTimes);
controlCorrectedApoptosis = zeros(nDoses,nTimes);
meanMTORCurves = zeros(nDoses,numel(tEval));
meanAutophagyCurves = zeros(nDoses,numel(tEval));

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% Sirolimus-only simulations
fprintf('============================================================\n');
fprintf('STEP 4 - SIROLIMUS-ONLY MECHANISM VALIDATION\n');
fprintf('Cisplatin C           : %.1f\n',C);
fprintf('Cells                 : %d\n',nCells);
fprintf('Random seed           : %d\n',randomSeed);
fprintf('kSiro assumption      : %.4f per uM\n',kSiro);
fprintf('Baseline apoptosis    : %d/%d (%.1f %%)\n', ...
    nBaseline,nCells,baselinePercent);
fprintf('============================================================\n');

for iDose = 1:nDoses
    Sirolimus = sirolimusDoses(iDose);
    deathTime = nan(nCells,1);
    mTORAll = zeros(nCells,numel(tEval));
    autophagyAll = zeros(nCells,numel(tEval));

    fprintf('Running sirolimus = %.1f uM (%d/%d) ...\n', ...
        Sirolimus,iDose,nDoses);

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

        autophagyAll(iCell,:) = yCell(:,2).';
        mTORAll(iCell,:) = yCell(:,7).';

        firstCrossing = find(yCell(:,10) >= BCLmit,1,'first');
        if ~isempty(firstCrossing)
            deathTime(iCell) = tCell(firstCrossing);
        end
    end

    meanMTORCurves(iDose,:) = mean(mTORAll,1);
    meanAutophagyCurves(iDose,:) = mean(autophagyAll,1);

    for iTime = 1:nTimes
        currentTime = evaluationTimes(iTime);
        timeIndex = find(abs(tEval-currentTime) < 1e-9,1,'first');

        meanMTORAtTime(iDose,iTime) = mean(mTORAll(:,timeIndex));
        meanAutophagyAtTime(iDose,iTime) = ...
            mean(autophagyAll(:,timeIndex));
        apoptosisAtTime(iDose,iTime) = ...
            100*mean(deathTime <= currentTime);
    end
end

%% Time-matched correction using sirolimus = 0 control
for iTime = 1:nTimes
    controlCorrectedApoptosis(:,iTime) = ...
        apoptosisAtTime(:,iTime)-apoptosisAtTime(1,iTime);
end

%% Results table
results = table( ...
    sirolimusDoses.', ...
    meanMTORAtTime(:,1),meanMTORAtTime(:,2), ...
    meanAutophagyAtTime(:,1),meanAutophagyAtTime(:,2), ...
    apoptosisAtTime(:,1),apoptosisAtTime(:,2), ...
    controlCorrectedApoptosis(:,1),controlCorrectedApoptosis(:,2), ...
    'VariableNames',{ ...
    'Sirolimus_uM', ...
    'MeanMTOR_24h','MeanMTOR_48h', ...
    'MeanAutophagy_24h','MeanAutophagy_48h', ...
    'Apoptosis_24h_pct','Apoptosis_48h_pct', ...
    'ControlCorrected_24h_pctPoint','ControlCorrected_48h_pctPoint'});

fprintf('\nRESULTS\n');
disp(results);

fprintf('Mechanism check, 0 to 20 uM:\n');
fprintf('  mTOR change at 48 h      : %.5f -> %.5f\n', ...
    meanMTORAtTime(1,2),meanMTORAtTime(end,2));
fprintf('  Autophagy change at 48 h : %.5f -> %.5f\n', ...
    meanAutophagyAtTime(1,2),meanAutophagyAtTime(end,2));

if meanMTORAtTime(end,2) < meanMTORAtTime(1,2) && ...
        meanAutophagyAtTime(end,2) > meanAutophagyAtTime(1,2)
    fprintf('Mechanism validation      : PASS\n');
else
    fprintf('Mechanism validation      : CHECK MODEL\n');
end

outputFile = 'Step4_SirolimusOnly_Results.csv';
writetable(results,outputFile);
fprintf('Results saved to: %s\n',outputFile);
fprintf('============================================================\n');

%% Figures
fig = figure('Color','w','Name','Step 4 - Sirolimus Mechanism Validation');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(sirolimusDoses,meanMTORAtTime(:,1),'-ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(sirolimusDoses,meanMTORAtTime(:,2),'-sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Sirolimus (uM)');
ylabel('Mean mTOR');
title('mTOR Inhibition');
legend('24 hours','48 hours','Location','best');
grid on;

nexttile;
plot(sirolimusDoses,meanAutophagyAtTime(:,1),'-ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(sirolimusDoses,meanAutophagyAtTime(:,2),'-sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Sirolimus (uM)');
ylabel('Mean ATPHG');
title('Autophagy Response');
legend('24 hours','48 hours','Location','best');
grid on;

nexttile;
plot(sirolimusDoses,apoptosisAtTime(:,1),'-ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(sirolimusDoses,apoptosisAtTime(:,2),'-sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Sirolimus (uM)');
ylabel('Apoptotic cells (%)');
title('Sirolimus-Only Apoptosis');
legend('24 hours','48 hours','Location','best');
ylim([0 100]);
grid on;

nexttile;
grayLevel = linspace(0,0.72,nDoses);
lineStyles = {'-','--',':','-.'};
hold on;
for iDose = 1:nDoses
    curveColor = grayLevel(iDose)*[1 1 1];
    styleIndex = mod(iDose-1,numel(lineStyles))+1;
    plot(tEval,meanAutophagyCurves(iDose,:), ...
        'Color',curveColor, ...
        'LineStyle',lineStyles{styleIndex}, ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('%.0f uM',sirolimusDoses(iDose)));
end
xlabel('Time (hours)');
ylabel('Mean ATPHG');
title('Autophagy Time Course');
legend('Location','eastoutside');
grid on;

% Explicit light theme.
allAxes = findall(fig,'Type','axes');
set(allAxes, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.70 0.70 0.70], ...
    'MinorGridColor',[0.85 0.85 0.85], ...
    'Box','on');
set(findall(fig,'Type','text'),'Color','k');

allLegends = findall(fig,'Type','legend');
set(allLegends,'Color','w','TextColor','k','EdgeColor','k');

end
