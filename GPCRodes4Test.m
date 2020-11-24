%% Kinetic Modeling of GPCR system 20.380

close all; clc; clear;

% parameters: 
V_polymersome = 0.5 * 10^-15; % in liters
N_av = 6.022 * 10^23; % avogadro's number

% initial concentration (all in M)
L = 1e-3                                       % odorant concentration
R = 10 / V_polymersome / N_av              % GPCR concentration
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
AC = 1.26*10^-5;                            % adenylyl cyclase
Ga_GTP_ACp = 0;                             % G alpha bound to GTP and activated AC

% cAMP = 10^-1;                               % cyclic AMP
% ATP = 3152*10^-6;                           % ATP
% Ga_GTP_ACp_ATP = 0;                         % G alpha bound to GTP, activated AC, and ATP

tspan = [0,10]; 

y0 = [L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp];

options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
[t,y] = ode23s(@GPCRodes4, tspan, y0, options);

%%

y = num2cell(y, 1);
[L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp] = y{:};

%%

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

%%

figure('Position', [10 10 1600 400]);
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
ylim([1e-12, R_tot(1)*10]);
hold off

%%

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
ylim([1e-12, Ga_tot(1)*10]);
hold off

%%

subplot(1,3,3);
semilogy(t,AC, 'LineWidth', 2)
hold on
semilogy(t,Ga_GTP_ACp, 'LineWidth', 2)
xlabel('Time (s)','FontSize',14)
ylabel('Concentration (M)', 'Fontsize', 14)
title(sprintf('AC Concentrations (Total %0.1e M)', AC_tot(1)),'Fontsize', 16)
legend('AC', 'G_\alpha-GTP:AC*');
ylim([1e-12, AC_tot(1)*10]);
hold off