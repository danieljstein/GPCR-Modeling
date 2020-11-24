%% Kinetic Modeling of Entire Polymersome system 20.380

close all; clc; clear;

%% parameters: 
V_polymersome = 0.5 * 10^-15; % in liters
N_av = 6.022 * 10^23; % avogadro's number

%% initial concentrations (all in M)

% GPCR
L = 1e-6                                       % odorant concentration
R = 10 / V_polymersome / N_av               % GPCR concentration
G = 100 / V_polymersome / N_av              % G with alpha, beta and gamma subunits bound with GDP
Rp = 0;                                     % GPCR active conformation
AR = 0;                                     % GPCR bound to ligand (inactive)
ARp = 0;                                    % GPCR bound to ligand (active)
RG = 0;                                     % GPCR bound to G protein (inactive)
RpG = 0;                                    % GPCR bound to G protein (active)
ARG = 0;                                    % GPCR bound to ligand and G protein (inactive)
ARpG = 0;                                   % GPCR bound to ligand and G protein (active)
Ga_GTP = 0;                                 % G alpha bound with GTP
Gby = 0;                                    % G beta and gamma
Ga_GDP = 0;                                 % G alpha bound to GTP
GTP = 468*10^-6;                            % GTP
GDP = 159*10^-6;                            % GDP
GTP_GDP_tot = G + GTP + GDP
AC = 1.26*10^-5                             % adenylyl cyclase
Ga_GTP_ACp = 0;                             % G alpha bound to GTP and activated AC
y_GPCR = [L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp];

% AC
E = Ga_GTP_ACp;                             % E (activated adenylyl cyclase - assumes unbound AC has no activity)
S = 3152*10^-6                              % ATP
ES = 0;                                     % E-ATP
EQP = 0;                                    % E-cAMP-PP
EQ = 0;                                     % E-cAMP
EP = 0;                                     % E-PP
Q = 0;                                      % cAMP
P = 0;                                      % PP
y_AC = [E, S, ES, EQP, EQ, EP, Q, P];

% PKA
A0B0 = 1.26*10^-4;                          % Protein Kinase A (PKA) regulatory subunit with catalytic subunit (assume [ATP] >> K_M_ATP = 4.22 µM)
                                            % PKA R:C with empty cAMP binding sites
A1B0 = 0;                                   % PKA R:C with cAMP A-site filled
A0B1 = 0;                                   % PKA R:C with cAMP B-site filled
A1B1 = 0;                                   % PKA R:C with both cAMP binding sites filled
cAMP = Q;                                   % cAMP
y_PKA = [A0B0, A1B0, A0B1, A1B1, cAMP];

S = 2.5e-3;                                 % peptide substrate for phosphorylation
P = 0;                                      % phosphorylated substrate
y_PKA_MM = [S, P];

%% ODEs

y0 = [y_GPCR y_AC y_PKA y_PKA_MM];
tspan = [0,60*5]; % time in seconds
options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
[t,y] = ode23s(@PolymersomeOdes, tspan, y0, options);

%%

y = num2cell(y, 1);
[L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp] = y{1:17};

% verify conservation of each component
L_tot = L + AR + ARp + ARG + ARpG;
L_tot = [mean(L_tot), std(L_tot)]

R_tot = R + Rp + AR + ARp + RG + RpG + ARG + ARpG;
R_tot = [mean(R_tot), std(R_tot)]

Ga_tot = G + RG + RpG + ARG + ARpG + Ga_GTP + Ga_GDP + Ga_GTP_ACp;
Ga_tot = [mean(Ga_tot), std(Ga_tot)]

Gby_tot = G + RG + RpG + ARG + ARpG + Gby;
Gby_tot = [mean(Gby_tot), std(Gby_tot)]

GDP_tot = G + RG + RpG + ARG + ARpG + Ga_GTP + Ga_GDP + GTP + GDP;
GDP_tot = [mean(GDP_tot), std(GDP_tot)]

AC_tot = AC + Ga_GTP_ACp;
AC_tot = [mean(AC_tot), std(AC_tot)]

figure('Position', [10 500 1600 400]);
subplot(1,3,1);
semilogy(t, R + Rp, 'LineWidth', 2)
hold on
semilogy(t, ARG + ARpG, 'LineWidth', 2)
semilogy(t, AR + ARp, 'LineWidth', 2)
semilogy(t, RG + RpG, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('GPCR Concentrations (Total %0.1e M)', R_tot(1)),'Fontsize', 16)
legend('GPCR', 'GPCR:L:G_{\alpha\beta\gamma}-GDP', 'GPCR:L', 'GPCR:G_{\alpha\beta\gamma}-GDP');
xlim([0,5]);
ylim([1e-12, R_tot(1)*10]);
hold off

subplot(1,3,2);
semilogy(t,G, 'LineWidth', 2)
hold on
semilogy(t,RG + RpG + ARG + ARpG, 'LineWidth', 2)
semilogy(t,Ga_GDP, 'LineWidth', 2)
semilogy(t,Ga_GTP, 'LineWidth', 2)
semilogy(t,Ga_GTP_ACp, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('G Protein Concentrations (Total %0.1e M)', Ga_tot(1)),'Fontsize', 16)
legend('G_{\alpha\beta\gamma}-GDP', 'GPCR*:G_{\alpha\beta\gamma}-GDP', 'G_\alpha-GDP', 'G_\alpha-GTP', 'G_\alpha-GTP:AC*')
xlim([0,5]);
ylim([1e-12, Ga_tot(1)*10]);
hold off

subplot(1,3,3);
semilogy(t,AC, 'LineWidth', 2)
hold on
semilogy(t,Ga_GTP_ACp, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('AC Concentrations (Total %0.1e M)', AC_tot(1)),'Fontsize', 16)
legend('AC', 'G_\alpha-GTP:AC*');
xlim([0,5]);
ylim([1e-12, AC_tot(1)*10]);
hold off

%%

[E, S, ES, EQP, EQ, EP, Q, P] = y{18:25};
[A0B0, A1B0, A0B1, A1B1, cAMP] = y{26:30};
[uPeptide, pPeptide] = y{31:32};

E_tot = E + ES + EQP + EQ + EP + AC;
E_tot = [mean(E_tot), std(E_tot)]
cAMP_tot = S + ES + EQ + Q + EQP + A1B0 + A0B1 + 2*A1B1;
cAMP_tot = [mean(cAMP_tot), std(cAMP_tot)]
P_tot = S + ES + EP + P + EQP;
P_tot = [mean(P_tot), std(P_tot)]

figure('Position', [10 10 800 400]);
subplot(1,2,1);
semilogy(t, E, 'LineWidth', 2)
hold on
semilogy(t, ES, 'LineWidth', 2)
semilogy(t, EQP, 'LineWidth', 2)
semilogy(t, EQ, 'LineWidth', 2)
semilogy(t, EP, 'LineWidth', 2)
semilogy(t, AC, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('AC Concentrations (Total %0.1e M)', E_tot(1)),'Fontsize', 16)
legend('AC*', 'AC*-ATP', 'AC*-cAMP-PP', 'AC*-cAMP', 'AC*-PP', 'AC');
xlim([0,5]);
ylim([1e-14, E_tot(1)*10]);
hold off

subplot(1,2,2);
semilogy(t, S, 'LineWidth', 2)
hold on
semilogy(t, ES, 'LineWidth', 2)
semilogy(t, EQP + EQ, 'LineWidth', 2)
semilogy(t, Q, 'LineWidth', 2)
semilogy(t, A1B0 + A0B1 + 2*A1B1, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('ATP/cAMP Concentrations (Total %0.1e M)', cAMP_tot(1)),'Fontsize', 16)
legend('ATP', 'AC*-ATP', 'AC*-cAMP-PP + AC*-cAMP', 'Free cAMP', 'cAMP bound to PKA');
ylim([1e-13, cAMP_tot(1)*10]);
hold off

%%

AB_tot = A0B0 + A1B0 + A0B1 + A1B1;
AB_tot = [mean(AB_tot), std(AB_tot)]

peptide_tot = uPeptide + pPeptide;
peptide_tot = [mean(peptide_tot), std(peptide_tot)];

figure('Position', [900 10 800 400]);
subplot(1,2,1)
semilogy(t, A0B0, 'DisplayName', 'A0B0', 'LineWidth', 2);
hold on;
semilogy(t, A1B1, 'DisplayName', 'A1B1', 'LineWidth', 2);
semilogy(t, A1B0, 'DisplayName', 'A1B0', 'LineWidth', 2);
semilogy(t, A0B1, 'DisplayName', 'A0B1', 'LineWidth', 2);
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('PKA Concentrations (Total %0.1e M)', AB_tot(1)),'Fontsize', 16)
legend();
ylim([1e-12, AB_tot(1)*10]);
hold off

subplot(1,2,2)
semilogy(t, uPeptide, 'DisplayName', 'Unphosphorylated', 'LineWidth', 2);
hold on;
semilogy(t, pPeptide, 'DisplayName', 'Phosphorylated', 'LineWidth', 2);
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('Peptide Substrate Concentrations (Total %0.1e M)', peptide_tot(1)),'Fontsize', 16)
legend();
ylim([1e-12, peptide_tot(1)*10]);
hold off