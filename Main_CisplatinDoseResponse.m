function Main_CisplatinDoseResponse
% STEP 3: Validate cisplatin dose-response in 100 virtual RPT cells.
% Put this file in the same folder as Fun.m, then run:
%   Main_CisplatinDoseResponse

% The same BCL2mit population is reused at every cisplatin dose so that
% differences between conditions are caused by dose, not random sampling.

clc;
close all;

if exist('Fun','file') ~= 2
    error(['Fun.m was not found. Put Main_CisplatinDoseResponse.m ' ...
           'in the same folder as Fun.m.']);
end

%% Simulation settings
cisplatinDoses = [0 2 4 5 6 10 20];
evaluationTimes = [24 48];             % hours
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

% No-drug steady-state initial conditions
y0 = [3 0.285 0.831 0.281 0.103 0.537 0.335 ...
      0.0184 0.397 0.069 0.0382 0.0221 0.717];

%% Generate one fixed population of 100 cells
% Target arithmetic mean and SD reported in Tavassoly et al.
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1 + (sdBCL2 / meanBCL2)^2));
muLog = log(meanBCL2) - 0.5 * sigmaLog^2;

rng(randomSeed,'twister');
BCLmitValues = exp(muLog + sigmaLog * randn(nCells,1));

baselineCells = BCLmitValues <= y0(10);
nBaseline = sum(baselineCells);
baselinePercent = 100 * nBaseline / nCells;

%% Allocate outputs
nDoses = numel(cisplatinDoses);
nTimes = numel(evaluationTimes);
rawApoptosis = zeros(nDoses,nTimes);
correctedApoptosis = zeros(nDoses,nTimes);
meanAutophagyAtTime = zeros(nDoses,nTimes);
medianDeathTime = nan(nDoses,1);
apoptosisCurves = zeros(nDoses,numel(tEval));

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% Dose-response simulation
fprintf('============================================================\n');
fprintf('STEP 3 - CISPLATIN DOSE-RESPONSE VALIDATION\n');
fprintf('Cells                 : %d\n',nCells);
fprintf('Random seed           : %d\n',randomSeed);
fprintf('Sample mean BCL2mit   : %.5f\n',mean(BCLmitValues));
fprintf('Sample SD BCL2mit     : %.5f\n',std(BCLmitValues));
fprintf('Baseline apoptosis    : %d/%d (%.1f %%)\n', ...
    nBaseline,nCells,baselinePercent);
fprintf('============================================================\n');

for iDose = 1:nDoses
    C = cisplatinDoses(iDose);
    deathTime = nan(nCells,1);
    autophagyAll = zeros(nCells,numel(tEval));

    fprintf('Running C = %.1f (%d/%d) ...\n',C,iDose,nDoses);

    for iCell = 1:nCells
        BCLmit = BCLmitValues(iCell);

        [tCell,yCell] = ode15s(@Fun,tEval,y0,options, ...
            ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl, ...
            sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal, ...
            kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0, ...
            wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m, ...
            wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP, ...
            Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp, ...
            sigmaA,wbh_c,wbh_a);

        autophagyAll(iCell,:) = yCell(:,2).';

        firstCrossing = find(yCell(:,10) >= BCLmit,1,'first');
        if ~isempty(firstCrossing)
            deathTime(iCell) = tCell(firstCrossing);
        end
    end

    % Cumulative apoptosis over the complete 48-hour simulation.
    for iPoint = 1:numel(tEval)
        apoptosisCurves(iDose,iPoint) = ...
            100 * mean(deathTime <= tEval(iPoint));
    end

    % Values at 24 and 48 hours.
    for iTime = 1:nTimes
        currentTime = evaluationTimes(iTime);
        timeIndex = find(abs(tEval-currentTime) < 1e-9,1,'first');

        rawApoptosis(iDose,iTime) = ...
            100 * mean(deathTime <= currentTime);
        meanAutophagyAtTime(iDose,iTime) = ...
            mean(autophagyAll(:,timeIndex));
    end

    % Median kematian hanya untuk sel yang awalnya masih hidup
    treatmentDeaths = deathTime( ...
        ~baselineCells & deathTime <= simulationTime);

    if ~isempty(treatmentDeaths)
        medianDeathTime(iDose) = median(treatmentDeaths);
    else
        medianDeathTime(iDose) = NaN;
    end
end

%% Correct every condition using the time-matched C = 0 control
for iTime = 1:nTimes
    correctedApoptosis(:,iTime) = ...
        max(0,rawApoptosis(:,iTime)-rawApoptosis(1,iTime));
end

%% Results table
results = table( ...
    cisplatinDoses.', ...
    rawApoptosis(:,1),rawApoptosis(:,2), ...
    correctedApoptosis(:,1),correctedApoptosis(:,2), ...
    meanAutophagyAtTime(:,1),meanAutophagyAtTime(:,2), ...
    medianDeathTime, ...
    'VariableNames',{ ...
    'Cisplatin_C', ...
    'Apoptosis_24h_pct','Apoptosis_48h_pct', ...
    'ControlCorrected_24h_pctPoint','ControlCorrected_48h_pctPoint', ...
    'MeanAutophagy_24h','MeanAutophagy_48h', ...
    'MedianDeathTime_h'});

fprintf('\nRESULTS\n');
disp(results);

outputFile = 'Step3_Cisplatin_DoseResponse_Results.csv';
writetable(results,outputFile);
fprintf('Results saved to: %s\n',outputFile);
fprintf('============================================================\n');

%% Figures
fig = figure('Color','w','Name','Step 3 - Cisplatin Dose-Response');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(cisplatinDoses,rawApoptosis(:,1),'-ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(cisplatinDoses,rawApoptosis(:,2),'-sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Cisplatin input, C');
ylabel('Apoptotic cells (%)');
title('Raw Dose-Response');
legend('24 hours','48 hours','Location','best');
ylim([0 100]);
grid on;

nexttile;
plot(cisplatinDoses,correctedApoptosis(:,1),'--ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(cisplatinDoses,correctedApoptosis(:,2),'--sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Cisplatin input, C');
ylabel('Increase over C = 0 (percentage points)');
title('Control-Corrected Apoptosis');
legend('24 hours','48 hours','Location','best');
ylim([0 100]);
grid on;

nexttile;
grayLevel = linspace(0,0.72,nDoses);
lineStyles = {'-','--',':','-.'};
hold on;
for iDose = 1:nDoses
    curveColor = grayLevel(iDose) * [1 1 1];
    styleIndex = mod(iDose-1,numel(lineStyles))+1;
    plot(tEval,apoptosisCurves(iDose,:), ...
        'Color',curveColor, ...
        'LineStyle',lineStyles{styleIndex}, ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('C = %.0f',cisplatinDoses(iDose)));
end
xlabel('Time (hours)');
ylabel('Apoptotic cells (%)');
title('Cumulative Apoptosis');
legend('Location','eastoutside');
ylim([0 100]);
grid on;

nexttile;
plot(cisplatinDoses,meanAutophagyAtTime(:,1),'-ok', ...
    'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6);
hold on;
plot(cisplatinDoses,meanAutophagyAtTime(:,2),'-sk', ...
    'LineWidth',2,'MarkerFaceColor','k','MarkerSize',6);
xlabel('Cisplatin input, C');
ylabel('Mean ATPHG');
title('Mean Autophagy');
legend('24 hours','48 hours','Location','best');
grid on;

% Explicit light theme: white background with black axes and text.
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
