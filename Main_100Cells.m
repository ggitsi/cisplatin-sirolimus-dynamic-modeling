function Main_100Cells
% STEP 2: Population simulation of 100 virtual RPT cells.
% Each cell uses one scalar BCL2mit value sampled from the lognormal
% distribution reported by Tavassoly et al. The same random seed should be
% reused in later treatment comparisons so that the virtual population is
% identical across conditions.

clc;
close all;

%% Simulation settings
C = 20;             % Cisplatin stress input; benchmark maps C=20 to 20 uM
simulationTime = 24;
nCells = 100;
randomSeed = 42;

%% Fixed model parameters from the original Main.m
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

% Original no-drug steady-state initial conditions
y0 = [3 0.285 0.831 0.281 0.103 0.537 0.335 ...
      0.0184 0.397 0.069 0.0382 0.0221 0.717];

%% Generate 100 BCL2mit values
% Arithmetic-space targets: mean=0.120 and SD=0.0292.
% lognrnd requires log-space parameters.
meanBCL2 = 0.120;
sdBCL2 = 0.0292;
sigmaLog = sqrt(log(1 + (sdBCL2 / meanBCL2)^2));
muLog = log(meanBCL2) - 0.5 * sigmaLog^2;

rng(randomSeed, 'twister');
% Equivalent to lognrnd(muLog,sigmaLog,[nCells,1]) but does not require
% the Statistics and Machine Learning Toolbox.
BCLmitValues = exp(muLog + sigmaLog * randn(nCells,1));

%% Run one ODE simulation per cell
tEval = linspace(0, simulationTime, 241); % output every 0.1 hour
options = odeset('RelTol',1e-3,'AbsTol',1e-3);

autophagyAll = zeros(nCells, numel(tEval));
bh3All = zeros(nCells, numel(tEval));
deathTime = nan(nCells,1);

for iCell = 1:nCells
    BCLmit = BCLmitValues(iCell); % scalar threshold for this cell

    [tCell,yCell] = ode15s(@Fun,tEval,y0,options, ...
        ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl, ...
        sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal, ...
        kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0, ...
        wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m, ...
        wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP, ...
        Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp, ...
        sigmaA,wbh_c,wbh_a);

    autophagyAll(iCell,:) = yCell(:,2).';
    bh3All(iCell,:) = yCell(:,10).';

    firstCrossing = find(yCell(:,10) >= BCLmit, 1, 'first');
    if ~isempty(firstCrossing)
        deathTime(iCell) = tCell(firstCrossing);
    end
end

apoptotic = ~isnan(deathTime);
nApoptotic = sum(apoptotic);
percentApoptosis = 100 * nApoptotic / nCells;
apoptosisCurve = zeros(size(tEval));

for iTime = 1:numel(tEval)
    apoptosisCurve(iTime) = 100 * mean(deathTime <= tEval(iTime));
end

meanAutophagy = mean(autophagyAll,1);
sdAutophagy = std(autophagyAll,0,1);

%% Command Window summary
fprintf('========================================\n');
fprintf('STEP 2 - 100 VIRTUAL CELLS\n');
fprintf('Cisplatin input C     : %.2f\n', C);
fprintf('Simulation time       : %.1f hours\n', simulationTime);
fprintf('Number of cells       : %d\n', nCells);
fprintf('Random seed           : %d\n', randomSeed);
fprintf('Sample mean BCL2mit   : %.5f\n', mean(BCLmitValues));
fprintf('Sample SD BCL2mit     : %.5f\n', std(BCLmitValues));
fprintf('Apoptotic cells       : %d/%d\n', nApoptotic, nCells);
fprintf('Final apoptosis       : %.1f %%\n', percentApoptosis);

if any(apoptotic)
    fprintf('Median death time     : %.2f hours\n', median(deathTime(apoptotic)));
    fprintf('Earliest death time   : %.2f hours\n', min(deathTime(apoptotic)));
    fprintf('Latest death time     : %.2f hours\n', max(deathTime(apoptotic)));
end
fprintf('========================================\n');

%% Population figures
fig = figure('Color','w','Name','Step 2 - Population of 100 Cells');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
histogram(BCLmitValues,12, ...
    'FaceColor',[0.35 0.35 0.35],'EdgeColor','k');
xline(meanBCL2,'--k','Target mean = 0.120','LineWidth',1.4);
xlabel('BCL2_{mit}');
ylabel('Number of cells');
title('Lognormal BCL2_{mit} Population');
grid on;

nexttile;
upperBand = meanAutophagy + sdAutophagy;
lowerBand = meanAutophagy - sdAutophagy;
fill([tEval fliplr(tEval)], [upperBand fliplr(lowerBand)], ...
    [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.70);
hold on;
plot(tEval,meanAutophagy,'k','LineWidth',2);
xlabel('Time (hours)');
ylabel('ATPHG');
title('Population Autophagy: Mean +/- SD');
legend('Mean +/- SD','Mean','Location','best');
grid on;

nexttile;
stairs(tEval,apoptosisCurve,'k','LineWidth',2);
xlabel('Time (hours)');
ylabel('Apoptotic cells (%)');
title('Cumulative Apoptosis');
ylim([0 100]);
grid on;

nexttile;
if any(apoptotic)
    histogram(deathTime(apoptotic),10, ...
        'FaceColor',[0.35 0.35 0.35],'EdgeColor','k');
    xlabel('Apoptosis time (hours)');
    ylabel('Number of cells');
    title('Distribution of Apoptosis Time');
else
    axis off;
    text(0.5,0.5,'No cells became apoptotic', ...
        'HorizontalAlignment','center','FontSize',12);
end
grid on;

% Force a light plot theme even when MATLAB itself uses dark mode.
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
set(allLegends, ...
    'Color','w', ...
    'TextColor','k', ...
    'EdgeColor','k');

end
