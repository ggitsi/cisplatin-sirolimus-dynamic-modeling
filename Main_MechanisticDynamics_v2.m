function Main_MechanisticDynamics_v2
% STEP 6: Mechanistic dynamics of cisplatin-sirolimus protection.
%
% Required file:
%   Fun_Sirolimus.m
%
% Main comparison:
%   Cisplatin C = 6
%   Sirolimus = 0 vs 10 uM
%   kSiro = 0.005 /uM (conservative nominal value used in Step 5)
%
% Important implementation detail:
%   Apoptosis is treated as an irreversible switch. Before the first
%   BH3 >= BCL2mit crossing, CASP = 0. At the first crossing, integration
%   is restarted with CASP = 1 and CASP remains 1 for the rest of the run.
%
% Outputs:
%   Step6_MechanisticDynamics_Endpoints.csv
%   Step6_MechanisticDynamics_Workspace.mat
%   Figure with mean mTOR, ATG13, ATPHG, Stress, BH3, and apoptosis.

clc;
close all;

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put Main_MechanisticDynamics.m ' ...
           'and Fun_Sirolimus.m in the same MATLAB folder.']);
end

%% Main Step 6 settings
C = 6;
sirolimusDoses = [0 10];      % uM
kSiro = 0.005;                % 1/uM, conservative nominal Step 5 value
evaluationTimes = [24 48];    % h
simulationTime = max(evaluationTimes);
timeStep = 0.1;
tEval = 0:timeStep:simulationTime;
nCells = 100;
randomSeed = 42;

%% Fixed Tavassoly parameters
p.BCL2T = 3;
p.IP3RT = 1;
p.CaT = 2;
p.wjnk_s = 1;
p.wdapk_s = 1;
p.wbcl_j = 1;
p.wmtor_ca = -1;
p.watg_m = -1;
p.wbh_a = 1;
p.wbec_d = 1;
p.Wcp_Ca = 1;
p.W_ATG5C = -1;
p.ka = 1.77;
p.kda = 0.0948;
p.ksb = 2.08;
p.kc = 0.510;
p.krb = 1.41;
p.kra = 3.83;
p.Gj = 1.73;
p.Gd = 4.05;
p.Gb = 5.21;
p.Gt = 1.61;
p.Gg = 1.04;
p.Gh = 0.0101;
p.Gl = 3.43;
p.Ga = 0.5240;
p.Gcp = 1.95;
p.sigmaj = 2.42;
p.sigmad = 1.01;
p.sigmab = 4.57;
p.sigmah = 2.89;
p.sigmag = 20.8;
p.sigmat = 3.51;
p.sigmal = 7.99;
p.sigmacp = 32.3;
p.sigmaA = 4.83;
p.kcasp = 2.01;
p.kin = 9.64;
p.kout = 6.31;
p.wjnk0 = -1.22;
p.wdapk0 = -2.98;
p.wbcl0 = -0.614;
p.wmtor0 = 0.202;
p.watg0 = 0.144;
p.wbh0 = -1.26;
p.wbec0 = -0.647;
p.W_ATG50 = 0.215;
p.Wcp_0 = -0.514;
p.wbh_s = 0.0801;
p.wbh_c = 0.0003;

% Original no-drug steady-state initial condition
p.y0 = [3 0.285 0.831 0.281 0.103 0.537 0.335 ...
        0.0184 0.397 0.069 0.0382 0.0221 0.717];

%% Same 100-cell BCL2mit population used in previous steps
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1 + (sdBCL2/meanBCL2)^2));
muLog = log(meanBCL2) - 0.5*sigmaLog^2;

rng(randomSeed,'twister');
BCLmitValues = exp(muLog + sigmaLog*randn(nCells,1));

%% Allocate population outputs
nCond = numel(sirolimusDoses);
nT = numel(tEval);

meanMTOR = zeros(nCond,nT);
meanATG13 = zeros(nCond,nT);
meanATPHG = zeros(nCond,nT);
meanStress = zeros(nCond,nT);
meanBH3 = zeros(nCond,nT);
apoptosisCurve = zeros(nCond,nT);
deathTimesAll = nan(nCells,nCond);

optionsBase = odeset('RelTol',1e-3,'AbsTol',1e-3);

fprintf('============================================================\n');
fprintf('STEP 6 - MECHANISTIC DYNAMICS\n');
fprintf('Cisplatin input C       : %.1f\n',C);
fprintf('Sirolimus comparison    : %.1f vs %.1f uM\n', ...
    sirolimusDoses(1),sirolimusDoses(2));
fprintf('kSiro                   : %.4f per uM\n',kSiro);
fprintf('Cells                   : %d\n',nCells);
fprintf('Random seed             : %d\n',randomSeed);
fprintf('Apoptosis switch        : irreversible\n');
fprintf('============================================================\n');

%% Run both treatment conditions
for iCond = 1:nCond
    Sirolimus = sirolimusDoses(iCond);

    mTORAll = zeros(nCells,nT);
    ATG13All = zeros(nCells,nT);
    ATPHGAll = zeros(nCells,nT);
    stressAll = zeros(nCells,nT);
    BH3All = zeros(nCells,nT);
    deathTime = nan(nCells,1);

    fprintf('Running C = %.1f, sirolimus = %.1f uM ...\n',C,Sirolimus);

    for iCell = 1:nCells
        BCLmit = BCLmitValues(iCell);

        [yCell,deathTime(iCell)] = simulateOneCellIrreversible( ...
            tEval,p,C,Sirolimus,kSiro,BCLmit,optionsBase);

        ATPHGAll(iCell,:) = yCell(:,2).';
        mTORAll(iCell,:) = yCell(:,7).';
        ATG13All(iCell,:) = yCell(:,8).';
        stressAll(iCell,:) = yCell(:,3).';
        BH3All(iCell,:) = yCell(:,10).';
    end

    deathTimesAll(:,iCond) = deathTime;
    meanMTOR(iCond,:) = mean(mTORAll,1);
    meanATG13(iCond,:) = mean(ATG13All,1);
    meanATPHG(iCond,:) = mean(ATPHGAll,1);
    meanStress(iCond,:) = mean(stressAll,1);
    meanBH3(iCond,:) = mean(BH3All,1);

    for iTime = 1:nT
        apoptosisCurve(iCond,iTime) = 100*mean(deathTime <= tEval(iTime));
    end
end

%% Endpoints at 24 and 48 hours
nTimes = numel(evaluationTimes);
rows = nCond*nTimes;
Condition = strings(rows,1);
Sirolimus_uM = zeros(rows,1);
Time_h = zeros(rows,1);
Mean_mTOR = zeros(rows,1);
Mean_ATG13 = zeros(rows,1);
Mean_ATPHG = zeros(rows,1);
Mean_Stress = zeros(rows,1);
Mean_BH3 = zeros(rows,1);
Apoptosis_pct = zeros(rows,1);

r = 0;
for iCond = 1:nCond
    for iTime = 1:nTimes
        r = r+1;
        idx = find(abs(tEval-evaluationTimes(iTime)) < 1e-12,1,'first');
        Condition(r) = sprintf('C=%.0f, Siro=%.0f uM',C,sirolimusDoses(iCond));
        Sirolimus_uM(r) = sirolimusDoses(iCond);
        Time_h(r) = evaluationTimes(iTime);
        Mean_mTOR(r) = meanMTOR(iCond,idx);
        Mean_ATG13(r) = meanATG13(iCond,idx);
        Mean_ATPHG(r) = meanATPHG(iCond,idx);
        Mean_Stress(r) = meanStress(iCond,idx);
        Mean_BH3(r) = meanBH3(iCond,idx);
        Apoptosis_pct(r) = apoptosisCurve(iCond,idx);
    end
end

endpoints = table(Condition,Sirolimus_uM,Time_h,Mean_mTOR,Mean_ATG13, ...
    Mean_ATPHG,Mean_Stress,Mean_BH3,Apoptosis_pct);

fprintf('\nENDPOINT SUMMARY\n');
disp(endpoints);

% Mechanistic change caused by sirolimus at 48 h
idx48 = find(abs(tEval-48) < 1e-12,1,'first');

fprintf('\nMECHANISTIC CHANGE AT 48 HOURS (10 uM vs 0 uM)\n');
fprintf('mTOR      : %.5f -> %.5f  (delta %.5f)\n', ...
    meanMTOR(1,idx48),meanMTOR(2,idx48), ...
    meanMTOR(2,idx48)-meanMTOR(1,idx48));
fprintf('ATG13     : %.5f -> %.5f  (delta %.5f)\n', ...
    meanATG13(1,idx48),meanATG13(2,idx48), ...
    meanATG13(2,idx48)-meanATG13(1,idx48));
fprintf('ATPHG     : %.5f -> %.5f  (delta %.5f)\n', ...
    meanATPHG(1,idx48),meanATPHG(2,idx48), ...
    meanATPHG(2,idx48)-meanATPHG(1,idx48));
fprintf('Stress    : %.5f -> %.5f  (delta %.5f)\n', ...
    meanStress(1,idx48),meanStress(2,idx48), ...
    meanStress(2,idx48)-meanStress(1,idx48));
fprintf('BH3       : %.5f -> %.5f  (delta %.5f)\n', ...
    meanBH3(1,idx48),meanBH3(2,idx48), ...
    meanBH3(2,idx48)-meanBH3(1,idx48));
fprintf('Apoptosis : %.1f%% -> %.1f%%  (protection %.1f percentage points)\n', ...
    apoptosisCurve(1,idx48),apoptosisCurve(2,idx48), ...
    apoptosisCurve(1,idx48)-apoptosisCurve(2,idx48));

%% Save outputs
outputFile = 'Step6_MechanisticDynamics_Endpoints.csv';
writetable(endpoints,outputFile);

save('Step6_MechanisticDynamics_Workspace.mat', ...
    'tEval','C','sirolimusDoses','kSiro','BCLmitValues','deathTimesAll', ...
    'meanMTOR','meanATG13','meanATPHG','meanStress','meanBH3', ...
    'apoptosisCurve','endpoints');

fprintf('\nResults saved to : %s\n',outputFile);
fprintf('Workspace saved  : Step6_MechanisticDynamics_Workspace.mat\n');
fprintf('============================================================\n');

%% Figure: mechanistic chain
fig = figure('Color','w','Name','Step 6 - Mechanistic Dynamics');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

labels = {'C=6, Siro=0 uM','C=6, Siro=10 uM'};

nexttile;
hold on;
plot(tEval,meanMTOR(1,:),'k-','LineWidth',2.0);
plot(tEval,meanMTOR(2,:),'k--','LineWidth',2.2);
xlabel('Time (hours)'); ylabel('Mean mTOR'); title('1. mTOR Inhibition');
legend(labels,'Location','best'); grid on;

nexttile;
hold on;
plot(tEval,meanATG13(1,:),'k-','LineWidth',2.0);
plot(tEval,meanATG13(2,:),'k--','LineWidth',2.2);
xlabel('Time (hours)'); ylabel('Mean ATG13'); title('2. ATG13 Activation');
legend(labels,'Location','best'); grid on;

nexttile;
hold on;
plot(tEval,meanATPHG(1,:),'k-','LineWidth',2.0);
plot(tEval,meanATPHG(2,:),'k--','LineWidth',2.2);
xlabel('Time (hours)'); ylabel('Mean ATPHG'); title('3. Autophagy');
legend(labels,'Location','best'); grid on;

nexttile;
hold on;
plot(tEval,meanStress(1,:),'k-','LineWidth',2.0);
plot(tEval,meanStress(2,:),'k--','LineWidth',2.2);
xlabel('Time (hours)'); ylabel('Mean stress, S'); title('4. Cellular Stress');
legend(labels,'Location','best'); grid on;

nexttile;
hold on;
plot(tEval,meanBH3(1,:),'k-','LineWidth',2.0);
plot(tEval,meanBH3(2,:),'k--','LineWidth',2.2);
yline(mean(BCLmitValues),':k','Mean BCL2_{mit}', ...
    'HandleVisibility','off','LineWidth',1.2);
xlabel('Time (hours)'); ylabel('Mean BH3'); title('5. Pro-apoptotic BH3');
legend(labels,'Location','best'); grid on;

nexttile;
hold on;
stairs(tEval,apoptosisCurve(1,:),'k-','LineWidth',2.0);
stairs(tEval,apoptosisCurve(2,:),'k--','LineWidth',2.2);
xlabel('Time (hours)'); ylabel('Apoptotic cells (%)'); title('6. Population Apoptosis');
legend(labels,'Location','best'); ylim([0 100]); grid on;

% Explicit light theme
allAxes = findall(fig,'Type','axes');
set(allAxes,'Color','w','XColor','k','YColor','k', ...
    'GridColor',[0.70 0.70 0.70], ...
    'MinorGridColor',[0.85 0.85 0.85], ...
    'Box','on');
set(findall(fig,'Type','text'),'Color','k');
allLegends = findall(fig,'Type','legend');
set(allLegends,'Color','w','TextColor','k','EdgeColor','k');

end

%% ------------------------------------------------------------------------
function [yAligned,deathTime] = simulateOneCellIrreversible( ...
    tEval,p,C,Sirolimus,kSiro,BCLmit,optionsBase)
% Simulate one cell with an irreversible apoptosis switch.

CASP_OFF = 0;
CASP_ON = 1;

y0 = p.y0;

% A small fraction of the generated population can already satisfy the
% threshold at t=0. Those cells start with CASP permanently ON.
if y0(10) >= BCLmit
    deathTime = 0;
    odeOn = @(t,y) modelRHS(t,y,p,C,Sirolimus,kSiro,BCLmit,CASP_ON);
    [~,yAligned] = ode15s(odeOn,tEval,y0,optionsBase);
    return;
end

% Phase 1: alive, CASP OFF, stop exactly at the first BH3-BCL2mit crossing.
odeOff = @(t,y) modelRHS(t,y,p,C,Sirolimus,kSiro,BCLmit,CASP_OFF);
eventFcn = @(t,y) apoptosisEvent(t,y,BCLmit);
optionsEvent = odeset(optionsBase,'Events',eventFcn);

[tPre,yPre,tEvent,yEvent] = ode15s(odeOff,tEval,y0,optionsEvent);

if isempty(tEvent)
    % No apoptosis during the simulated interval.
    deathTime = NaN;
    yAligned = yPre;
    return;
end

% First threshold crossing = irreversible commitment to apoptosis.
deathTime = tEvent(1);
ySwitch = yEvent(1,:).';

% Phase 2: CASP permanently ON after commitment.
tRemaining = tEval(tEval > deathTime + 1e-10);

if isempty(tRemaining)
    tAll = tPre;
    yAll = yPre;
else
    odeOn = @(t,y) modelRHS(t,y,p,C,Sirolimus,kSiro,BCLmit,CASP_ON);
    tSpanPost = [deathTime tRemaining];
    [tPost,yPost] = ode15s(odeOn,tSpanPost,ySwitch,optionsBase);

    % tPre ends at the event and tPost starts at the same event.
    tAll = [tPre; tPost(2:end)];
    yAll = [yPre; yPost(2:end,:)];
end

% All requested grid points are present except that the event introduces one
% extra off-grid time point. Interpolation simply removes that extra point
% and returns values on the common population time grid.
yAligned = interp1(tAll,yAll,tEval(:),'linear');

end

%% ------------------------------------------------------------------------
function dy = modelRHS(t,y,p,C,Sirolimus,kSiro,BCLmit,CASP)
% Wrapper around the user's Fun_Sirolimus.m using a compact parameter struct.

dy = Fun_Sirolimus(t,y, ...
    p.ka,p.kda,C,p.kc,p.ksb,p.krb,p.kra, ...
    p.Gj,p.Gd,p.Gb,p.Gt,p.Gg,p.Gh,p.Gl, ...
    p.sigmaj,p.sigmad,p.sigmab,p.sigmat,p.sigmag,p.sigmah,p.sigmal, ...
    p.kcasp,p.kin,p.kout,p.wbec0,p.wjnk0,p.wdapk0,p.wbcl0,p.watg0, ...
    p.wmtor0,p.wbh0,p.wjnk_s,p.wdapk_s,p.wbec_d,p.wbcl_j,p.watg_m, ...
    p.wmtor_ca,p.wbh_s,p.BCL2T,p.IP3RT,BCLmit,p.CaT,CASP, ...
    p.Wcp_0,p.Wcp_Ca,p.W_ATG50,p.W_ATG5C,p.Ga,p.Gcp,p.sigmacp, ...
    p.sigmaA,p.wbh_c,p.wbh_a,Sirolimus,kSiro);

end

%% ------------------------------------------------------------------------
function [value,isterminal,direction] = apoptosisEvent(~,y,BCLmit)
% Detect the first upward crossing of BH3 over the cell-specific BCL2mit.
value = y(10)-BCLmit;
isterminal = 1;
direction = 1;
end
