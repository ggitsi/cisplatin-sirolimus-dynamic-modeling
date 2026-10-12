
function Main
C=20; BCL2T=3; IP3RT=1; CaT=2; CASP =0;  
wjnk_s = 1; wdapk_s = 1; wbcl_j = 1; wmtor_ca = -1; watg_m = -1; wbh_a = 1; wbec_d = 1; Wcp_Ca = 1; W_ATG5C = -1; ka = 1.77; kda = 0.0948; ksb = 2.08; kc = 0.510; krb = 1.41; kra = 3.83; Gj = 1.73; Gd = 4.05; Gb = 5.21; Gt = 1.61; Gg = 1.04; Gh = 0.0101; Gl = 3.43; Ga = 0.5240; Gcp = 1.95; sigmaj = 2.42; sigmad = 1.01; sigmab = 4.57; sigmah = 2.89; sigmag = 20.8; sigmat = 3.51; sigmal = 7.99; sigmacp = 32.3; sigmaA = 4.83; kcasp = 2.01; kin = 9.64; 
kout = 6.31; wjnk0 = -1.22; wdapk0 = -2.98; wbcl0 = -0.614; wmtor0 = 0.202; watg0 = 0.144; wbh0 = -1.26; wbec0 = -0.647; W_ATG50 = 0.215; Wcp_0 = -0.514; wbh_s = 0.0801; wbh_c = 0.0003; 
 
y0 = [3 0.285 0.831 0.281 0.103 0.537 0.335 0.0184 0.397 0.069 0.0382 0.0221 0.717]; % Initial steady state values


BCLmit = 0.12; % lognormal disribution of BCL2mit

options = odeset('RelTol',1e-3,'AbsTol',1e-3); 
[t,y] = ode15s(@Fun,[0 24],y0,options,ka,kda,C,kc,ksb,krb,kra,Gj,Gd,Gb,Gt,Gg,Gh,Gl,sigmaj,sigmad,sigmab,sigmat,sigmag,sigmah,sigmal,kcasp,kin,kout,wbec0,wjnk0,wdapk0,wbcl0,watg0,wmtor0,wbh0,wjnk_s,wdapk_s,wbec_d,wbcl_j,watg_m,wmtor_ca,wbh_s,BCL2T,IP3RT,BCLmit,CaT,CASP,Wcp_0,Wcp_Ca,W_ATG50,W_ATG5C,Ga,Gcp,sigmacp,sigmaA,wbh_c,wbh_a);
clc;
close all;

% Membuat indikator CASP berdasarkan crossing BH3-BCLmit
CASP_plot = cummax(double(y(:,10) >= BCLmit));

% Mencari waktu pertama kali apoptosis terjadi
idxApoptosis = find(CASP_plot == 1, 1, 'first');

fprintf('========================================\n');
fprintf('BASELINE SINGLE-CELL SIMULATION\n');
fprintf('Cisplatin input C : %.2f\n', C);
fprintf('BCL2mit           : %.3f\n', BCLmit);
fprintf('Simulation time   : %.1f hours\n', t(end));

if isempty(idxApoptosis)
    fprintf('Result            : No apoptosis\n');
else
    fprintf('Result            : Apoptosis occurred\n');
    fprintf('Apoptosis time    : %.2f hours\n', t(idxApoptosis));
end

fprintf('========================================\n');

figure('Color','w','Name','Baseline Cisplatin Simulation');

tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% Grafik autophagy
nexttile;
plot(t,y(:,2),'b','LineWidth',2);
xlabel('Time (hours)');
ylabel('ATPHG');
title('Autophagy');
grid on;

% Grafik stres
nexttile;
plot(t,y(:,3),'Color',[0.85 0.33 0.10],'LineWidth',2);
xlabel('Time (hours)');
ylabel('Stress');
title('Cellular Stress');
grid on;

% Grafik BH3 dan BCL2mit
nexttile;
plot(t,y(:,10),'m','LineWidth',2);
hold on;
yline(BCLmit,'--k','BCL2_{mit}','LineWidth',1.5);
xlabel('Time (hours)');
ylabel('Concentration');
title('Apoptosis Threshold');
legend('BH3','BCL2_{mit}','Location','best');
grid on;

% Grafik CASP
nexttile;
stairs(t,CASP_plot,'r','LineWidth',2);
xlabel('Time (hours)');
ylabel('CASP');
title('Apoptosis State');
ylim([-0.05 1.05]);
yticks([0 1]);
yticklabels({'Alive','Apoptosis'});
grid on;


    