close all;
clc;
clear;

% Test AC ODEs
%
% Model from https://www.jbc.org/content/272/44/27787.full.pdf
% Recapitulates experimental results from Figure 1

velinv = [0, 0, 0, 0];
Sinv = [0, 0, 0, 0];
allS = [0.03e-3, 0.1e-3, 0.3e-3, 1e-3];


for jj = 1:4

    tspan = [0,1000]; 

    E = 0.7e-9;   % E
    S = allS(jj);   % ATP
    ES = 0;  % E-ATP
    EQP = 0; % E-cAMP-PP
    EQ = 0;  % E-cAMP
    EP = 0;  % E-PP
    Q = 0;   % cAMP
    P = 0;   % PP

    y0 = [E, S, ES, EQP, EQ, EP, Q, P];

    options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
    [t,y] = ode23s(@ACodes, tspan, y0, options);

    E = y(:,1);   % E
    S = y(:,2);   % ATP
    ES = y(:,3);  % E-ATP
    EQP = y(:,4); % E-cAMP-PP
    EQ = y(:,5);  % E-cAMP
    EP = y(:,6);  % E-PP
    Q = y(:,7);   % cAMP
    P = y(:,8);   % PP

    E0 = E + ES + EQP + EQ + EP;
    E0 = [mean(E0), std(E0)]
    Q0 = S + ES + EQ + Q + EQP;
    Q0 = [mean(Q0), std(Q0)]
    P0 = S + ES + EP + P + EQP;
    P0 = [mean(P0), std(P0)]

    name = {'E', 'S', 'ES', 'EQP', 'EQ', 'EP', 'Q', 'P'};
    
    figure;
    for ii = 1:8
        subplot(2,4,ii);
        plot(t, y(:,ii), 'LineWidth', 2)
        xlabel('Time (s)','FontSize',14)
        ylabel(name(ii), 'Fontsize', 14)
        title(allS(jj));
    end

    vel = (Q(end) - Q(1))*1e6 / (tspan(2)/60) / (mean(E0) * 60000 * 1000); % vel (µmol/min/mg), guessing molar mass of VC1 domain of AC
    velinv(jj) = 1/vel;
    Sinv(jj) = 1/(allS(jj)*1e6);
end

figure;
plot(Sinv, velinv, 'LineWidth', 2);
xlabel('1/[ATP] (µM)^{-1}');
ylabel('1/Velocity (µmol/min-mg)^{-1}');
% xlim([0,0.04]);
ylim([0,0.2]);
title('Simulated Adenylyl Cyclase Activity');