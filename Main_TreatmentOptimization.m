function Main_TreatmentOptimization
% ============================================================
% STEP 7 - NEPHROTOXICITY-CONSTRAINED TREATMENT OPTIMIZATION
%
% Goal:
% Evaluate cisplatin-sirolimus combinations at 48 h and determine
% the highest cisplatin model input C that remains below selected
% apoptosis thresholds.
%
% IMPORTANT:
% C is the MODEL INPUT used in the Tavassoly model.
% It must NOT automatically be interpreted as cisplatin dose in uM.
%
% Required file:
%   Fun_Sirolimus.m
%
% Outputs:
%   Step7_ApoptosisGrid_48h.csv
%   Step7_NephrotoxicityOptimization_48h.csv
%   Step7_MinSirolimus_Target20pct.csv
%   Step7_TreatmentOptimization_Workspace.mat
% ============================================================

clc;
close all;

%% ============================================================
% CHECK REQUIRED FUNCTION
% =============================================================

if exist('Fun_Sirolimus','file') ~= 2
    error(['Fun_Sirolimus.m was not found. Put this file and ' ...
           'Fun_Sirolimus.m in the same MATLAB folder.']);
end

%% ============================================================
% OUTPUT FOLDER
% =============================================================

outputFolder = ...
    'D:\Kuliah\Transfer Kredit\Final Project Dynamic Modelling';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% ============================================================
% STEP 7 SETTINGS
% =============================================================

% Cisplatin MODEL INPUT grid
cisplatinGrid = 0:0.5:20;

% Sirolimus concentrations
sirolimusDoses = [0 2 5 10 20];     % uM

% Conservative nominal scenario from Step 5
kSiro = 0.005;                       % per uM

simulationTime = 48;                 % hours

nCells = 100;
randomSeed = 42;

% User-defined nephrotoxicity thresholds
targetApoptosis = [10 20 30];        % %

% Primary threshold used for minimum-sirolimus analysis
primaryTarget = 20;                  % %

%% ============================================================
% FIXED PARAMETERS FROM TAVASSOLY MODEL
% =============================================================

BCL2T = 3;
IP3RT = 1;
CaT = 2;
CASP = 0;

wjnk_s = 1;
wdapk_s = 1;
wbcl_j = 1;
wmtor_ca = -1;
watg_m = -1;
wbh_a = 1;
wbec_d = 1;
Wcp_Ca = 1;
W_ATG5C = -1;

ka = 1.77;
kda = 0.0948;
ksb = 2.08;
kc = 0.510;
krb = 1.41;
kra = 3.83;

Gj = 1.73;
Gd = 4.05;
Gb = 5.21;
Gt = 1.61;
Gg = 1.04;
Gh = 0.0101;
Gl = 3.43;
Ga = 0.5240;
Gcp = 1.95;

sigmaj = 2.42;
sigmad = 1.01;
sigmab = 4.57;
sigmah = 2.89;
sigmag = 20.8;
sigmat = 3.51;
sigmal = 7.99;
sigmacp = 32.3;
sigmaA = 4.83;

kcasp = 2.01;
kin = 9.64;
kout = 6.31;

wjnk0 = -1.22;
wdapk0 = -2.98;
wbcl0 = -0.614;
wmtor0 = 0.202;
watg0 = 0.144;
wbh0 = -1.26;
wbec0 = -0.647;

W_ATG50 = 0.215;
Wcp_0 = -0.514;

wbh_s = 0.0801;
wbh_c = 0.0003;

%% ============================================================
% INITIAL CONDITIONS
% =============================================================

y0 = [ ...
    3 ...       % BECNT
    0.285 ...   % ATPHG
    0.831 ...   % Stress
    0.281 ...   % JNK
    0.103 ...   % DAPK
    0.537 ...   % BCLP
    0.335 ...   % mTOR
    0.0184 ...  % ATG13
    0.397 ...   % Ca
    0.069 ...   % BH3
    0.0382 ...  % BECNP
    0.0221 ...  % calpain
    0.717];     % ATG5

%% ============================================================
% GENERATE SAME 100-CELL POPULATION
% =============================================================

meanBCL2 = 0.120;
sdBCL2 = 0.0292;

sigmaLog = sqrt(log(1 + (sdBCL2 / meanBCL2)^2));

muLog = log(meanBCL2) ...
    - 0.5 * sigmaLog^2;

rng(randomSeed,'twister');

BCLmitValues = ...
    exp(muLog + sigmaLog * randn(nCells,1));

baselineCells = BCLmitValues <= y0(10);

baselineApoptosis = ...
    100 * mean(baselineCells);

%% ============================================================
% ALLOCATE RESULTS
% =============================================================

nC = numel(cisplatinGrid);
nS = numel(sirolimusDoses);

apoptosis = zeros(nC,nS);

optionsBase = ...
    odeset('RelTol',1e-3,'AbsTol',1e-3);

totalConditions = nC * nS;

conditionCounter = 0;

fprintf('============================================================\n');
fprintf('STEP 7 - NEPHROTOXICITY-CONSTRAINED OPTIMIZATION\n');
fprintf('Cells                    : %d\n',nCells);
fprintf('Random seed              : %d\n',randomSeed);
fprintf('Simulation time          : %.1f h\n',simulationTime);
fprintf('kSiro                    : %.4f per uM\n',kSiro);
fprintf('Cisplatin C range        : %.1f to %.1f\n', ...
    min(cisplatinGrid),max(cisplatinGrid));
fprintf('Cisplatin grid step      : %.2f\n', ...
    cisplatinGrid(2)-cisplatinGrid(1));
fprintf('Sirolimus doses          :');
fprintf(' %.1f',sirolimusDoses);
fprintf(' uM\n');
fprintf('Baseline apoptosis       : %.1f %%\n', ...
    baselineApoptosis);

fprintf('Target apoptosis levels  :');
fprintf(' %.0f%%',targetApoptosis);
fprintf('\n');

fprintf('Total conditions         : %d\n', ...
    totalConditions);

fprintf('Total maximum ODE runs   : %d\n', ...
    totalConditions*nCells);

fprintf('============================================================\n');

%% ============================================================
% RUN OPTIMIZATION GRID
% =============================================================

for iC = 1:nC

    C = cisplatinGrid(iC);

    for iS = 1:nS

        Sirolimus = sirolimusDoses(iS);

        conditionCounter = conditionCounter + 1;

        fprintf(['Condition %d/%d: C = %.2f, ' ...
                 'Sirolimus = %.1f uM\n'], ...
                 conditionCounter,totalConditions, ...
                 C,Sirolimus);

        dead = false(nCells,1);

        for iCell = 1:nCells

            BCLmit = BCLmitValues(iCell);

            % Cell already above threshold at baseline
            if y0(10) >= BCLmit
                dead(iCell) = true;
                continue;
            end

            % Stop integration when BH3 reaches BCL2mit.
            eventFunction = ...
                @(t,y) apoptosisEvent(t,y,BCLmit);

            optionsCell = ...
                odeset(optionsBase, ...
                'Events',eventFunction);

            odeFunction = @(t,y) ...
                Fun_Sirolimus( ...
                t,y, ...
                ka,kda,C,kc,ksb,krb,kra, ...
                Gj,Gd,Gb,Gt,Gg,Gh,Gl, ...
                sigmaj,sigmad,sigmab,sigmat, ...
                sigmag,sigmah,sigmal, ...
                kcasp,kin,kout, ...
                wbec0,wjnk0,wdapk0,wbcl0, ...
                watg0,wmtor0,wbh0, ...
                wjnk_s,wdapk_s,wbec_d, ...
                wbcl_j,watg_m,wmtor_ca,wbh_s, ...
                BCL2T,IP3RT,BCLmit,CaT,CASP, ...
                Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C, ...
                Ga,Gcp,sigmacp,sigmaA, ...
                wbh_c,wbh_a, ...
                Sirolimus,kSiro);

            [~,~,te] = ...
                ode15s( ...
                odeFunction, ...
                [0 simulationTime], ...
                y0, ...
                optionsCell);

            if ~isempty(te)
                dead(iCell) = true;
            end

        end

        apoptosis(iC,iS) = ...
            100 * mean(dead);

    end
end

%% ============================================================
% CALCULATE PROTECTION
% =============================================================

% Cisplatin-only apoptosis is the first sirolimus column.
monoApoptosis = apoptosis(:,1);

protection = zeros(size(apoptosis));

for iS = 1:nS
    protection(:,iS) = ...
        monoApoptosis - apoptosis(:,iS);
end

%% ============================================================
% FULL GRID TABLE
% =============================================================

[cGrid,sGrid] = ...
    ndgrid(cisplatinGrid,sirolimusDoses);

gridTable = table( ...
    cGrid(:), ...
    sGrid(:), ...
    repmat(kSiro,numel(cGrid),1), ...
    apoptosis(:), ...
    protection(:), ...
    'VariableNames',{ ...
    'Cisplatin_C', ...
    'Sirolimus_uM', ...
    'kSiro_per_uM', ...
    'Apoptosis_48h_pct', ...
    'Protection_vs_CisplatinOnly_pctPoint'});

gridFile = fullfile( ...
    outputFolder, ...
    'Step7_ApoptosisGrid_48h.csv');

writetable(gridTable,gridFile);

%% ============================================================
% MAXIMUM C UNDER EACH APOPTOSIS THRESHOLD
% =============================================================

nTargets = numel(targetApoptosis);

nRows = nTargets*nS;

thresholdOut = zeros(nRows,1);
siroOut = zeros(nRows,1);
maxCOut = nan(nRows,1);
gainOut = nan(nRows,1);
upperBoundOut = false(nRows,1);

iRow = 0;

for iTarget = 1:nTargets

    target = targetApoptosis(iTarget);

    % Maximum C without sirolimus
    idxMono = find( ...
        apoptosis(:,1) <= target, ...
        1,'last');

    if isempty(idxMono)
        maxCMono = NaN;
    else
        maxCMono = cisplatinGrid(idxMono);
    end

    for iS = 1:nS

        iRow = iRow + 1;

        idx = find( ...
            apoptosis(:,iS) <= target, ...
            1,'last');

        thresholdOut(iRow) = target;
        siroOut(iRow) = sirolimusDoses(iS);

        if ~isempty(idx)

            maxCOut(iRow) = ...
                cisplatinGrid(idx);

            gainOut(iRow) = ...
                maxCOut(iRow) - maxCMono;

            % If the last grid point is still acceptable,
            % true maximum may be outside scanned range.
            if idx == nC
                upperBoundOut(iRow) = true;
            end

        end
    end
end

optimizationTable = table( ...
    thresholdOut, ...
    siroOut, ...
    maxCOut, ...
    gainOut, ...
    upperBoundOut, ...
    'VariableNames',{ ...
    'ApoptosisThreshold_pct', ...
    'Sirolimus_uM', ...
    'MaxC_withinScannedGrid', ...
    'IncreaseInToleratedC_vs_NoSirolimus', ...
    'ReachedScanUpperBound'});

optimizationFile = fullfile( ...
    outputFolder, ...
    'Step7_NephrotoxicityOptimization_48h.csv');

writetable(optimizationTable,optimizationFile);

%% ============================================================
% MINIMUM SIROLIMUS REQUIRED FOR PRIMARY TARGET
% =============================================================

minSiro = nan(nC,1);
targetReached = false(nC,1);

for iC = 1:nC

    idx = find( ...
        apoptosis(iC,:) <= primaryTarget, ...
        1,'first');

    if ~isempty(idx)

        minSiro(iC) = ...
            sirolimusDoses(idx);

        targetReached(iC) = true;

    end
end

minimumSiroTable = table( ...
    cisplatinGrid.', ...
    minSiro, ...
    targetReached, ...
    'VariableNames',{ ...
    'Cisplatin_C', ...
    'MinimumSirolimus_uM', ...
    'TargetAchieved'});

minimumSiroFile = fullfile( ...
    outputFolder, ...
    sprintf( ...
    'Step7_MinSirolimus_Target%dpct.csv', ...
    primaryTarget));

writetable(minimumSiroTable,minimumSiroFile);

%% ============================================================
% COMMAND WINDOW RESULTS
% =============================================================

fprintf('\n============================================================\n');
fprintf('STEP 7 OPTIMIZATION SUMMARY\n');
fprintf('============================================================\n');

disp(optimizationTable);

fprintf('\nPRIMARY TARGET = %.0f%% APOPTOSIS\n', ...
    primaryTarget);

primaryRows = ...
    optimizationTable.ApoptosisThreshold_pct ...
    == primaryTarget;

disp(optimizationTable(primaryRows,:));

fprintf('\nInterpretation:\n');
fprintf(['MaxC_withinScannedGrid = highest cisplatin MODEL INPUT C ' ...
         'with apoptosis <= selected threshold.\n']);

fprintf(['IncreaseInToleratedC_vs_NoSirolimus = shift in the ' ...
         'nephrotoxicity-constrained C window.\n']);

fprintf(['This is NOT proof of improved antitumor efficacy and ' ...
         'C must not automatically be interpreted as uM cisplatin.\n']);

%% ============================================================
% SAVE WORKSPACE
% =============================================================

workspaceFile = fullfile( ...
    outputFolder, ...
    'Step7_TreatmentOptimization_Workspace.mat');

save(workspaceFile, ...
    'cisplatinGrid', ...
    'sirolimusDoses', ...
    'kSiro', ...
    'simulationTime', ...
    'targetApoptosis', ...
    'primaryTarget', ...
    'BCLmitValues', ...
    'baselineApoptosis', ...
    'apoptosis', ...
    'protection', ...
    'gridTable', ...
    'optimizationTable', ...
    'minimumSiroTable');

%% ============================================================
% FIGURES
% =============================================================

fig = figure( ...
    'Color','w', ...
    'Name','Step 7 - Treatment Optimization');

tiledlayout(2,2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

%% ------------------------------------------------------------
% PANEL 1: APOPTOSIS VS CISPLATIN C
% ------------------------------------------------------------

nexttile;
hold on;

for iS = 1:nS

    plot( ...
        cisplatinGrid, ...
        apoptosis(:,iS), ...
        'LineWidth',1.8, ...
        'Marker','o', ...
        'MarkerSize',4, ...
        'DisplayName', ...
        sprintf('Siro = %.0f uM', ...
        sirolimusDoses(iS)));

end

yline( ...
    primaryTarget, ...
    '--', ...
    sprintf('Target = %.0f%%',primaryTarget), ...
    'HandleVisibility','off');

xlabel('Cisplatin model input, C');
ylabel('Apoptotic cells at 48 h (%)');

title('48-h Nephrotoxicity Response');

ylim([0 100]);
grid on;
legend('Location','best');

%% ------------------------------------------------------------
% PANEL 2: HEATMAP
% ------------------------------------------------------------

nexttile;

imagesc( ...
    sirolimusDoses, ...
    cisplatinGrid, ...
    apoptosis);

set(gca,'YDir','normal');

xlabel('Sirolimus (uM)');
ylabel('Cisplatin model input, C');

title('Apoptosis at 48 h (%)');

colorbar;
caxis([0 100]);

%% ------------------------------------------------------------
% PANEL 3: MAXIMUM C UNDER EACH THRESHOLD
% ------------------------------------------------------------

nexttile;
hold on;

for iTarget = 1:nTargets

    target = targetApoptosis(iTarget);

    rows = ...
        optimizationTable.ApoptosisThreshold_pct ...
        == target;

    plot( ...
        optimizationTable.Sirolimus_uM(rows), ...
        optimizationTable.MaxC_withinScannedGrid(rows), ...
        'LineWidth',1.8, ...
        'Marker','o', ...
        'MarkerSize',6, ...
        'DisplayName', ...
        sprintf('Apoptosis <= %.0f%%',target));

end

xlabel('Sirolimus (uM)');
ylabel('Maximum cisplatin model input, C');

title('Nephrotoxicity-Constrained C Window');

grid on;
legend('Location','best');

%% ------------------------------------------------------------
% PANEL 4: PROTECTION AT KEY C VALUES
% ------------------------------------------------------------

nexttile;
hold on;

focusC = [6 10 20];

for iFocus = 1:numel(focusC)

    idxC = find( ...
        abs(cisplatinGrid-focusC(iFocus)) < 1e-9, ...
        1,'first');

    if ~isempty(idxC)

        plot( ...
            sirolimusDoses, ...
            protection(idxC,:), ...
            'LineWidth',1.8, ...
            'Marker','o', ...
            'MarkerSize',6, ...
            'DisplayName', ...
            sprintf('C = %.0f',focusC(iFocus)));

    end
end

yline(0,'--','HandleVisibility','off');

xlabel('Sirolimus (uM)');
ylabel('Protection at 48 h (percentage points)');

title('Protection vs Cisplatin Alone');

grid on;
legend('Location','best');

%% ============================================================
% LIGHT FIGURE THEME
% =============================================================

allAxes = findall(fig,'Type','axes');

set(allAxes, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.75 0.75 0.75], ...
    'Box','on');

set(findall(fig,'Type','text'),'Color','k');

%% ============================================================
% SAVE FIGURE
% =============================================================

figureFile = fullfile( ...
    outputFolder, ...
    'Step7_TreatmentOptimization_Figure.png');

try
    exportgraphics( ...
        fig, ...
        figureFile, ...
        'Resolution',300);
catch
    fprintf('Figure could not be exported automatically.\n');
end

%% ============================================================
% FINAL OUTPUT PATHS
% =============================================================

fprintf('\n============================================================\n');
fprintf('FILES SAVED\n');
fprintf('============================================================\n');

fprintf('Grid results       : %s\n',gridFile);
fprintf('Optimization       : %s\n',optimizationFile);
fprintf('Minimum sirolimus  : %s\n',minimumSiroFile);
fprintf('Workspace          : %s\n',workspaceFile);
fprintf('Figure             : %s\n',figureFile);

fprintf('============================================================\n');

end


%% ============================================================
% LOCAL EVENT FUNCTION
% =============================================================

function [value,isterminal,direction] = ...
    apoptosisEvent(~,y,BCLmit)

% Apoptosis commitment occurs when:
% BH3 >= BCL2mit

value = y(10) - BCLmit;

% Stop integration at the first crossing
isterminal = 1;

% We only need an upward crossing
direction = 1;

end