close all;
clc;
clear;

% Test PKA Michaelis Menten ODEs
%
% Rate constants from http://www.perkinelmer.co.jp/Portals/0/resource/products_ls/reader/pdf/LC3000-AP-207.pdf
% Recapitulates experimental results from Figures 2 and 4
% Assumes [ATP] >> K_M_ATP = 4.22 µM (see Figure 3)

E0 = [2.0, 0.67, 0.22, 0.074, 0.025]*1e-9;  % M
N = length(E0);

figure(1);
hold on;

for jj = 1:N

    tspan = [0,3600];

    S = 1e-6;   % M
    P = 0;      % M

    y0 = [S, P];

    options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
    [t,y] = ode23s(@PKAMichaelisMentenOdes, tspan, y0, options, E0(jj));
    
    S = y(:,1);
    P = y(:,2);
    
%     S_P_tot = S + P;    % M
%     mean(S_P_tot)
%     var(S_P_tot)
    
    plot(t/60, P ./ (P + S), 'LineWidth', 2, 'DisplayName', num2str(E0(jj)));
end

ylabel('Time (min)');
xlabel('P/P+S');
legend();