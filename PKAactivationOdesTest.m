close all;
clc;
clear;

% Test PKA activation ODEs
%
% Binding constants from https://www.nature.com/articles/s41467-019-11930-2
% Recapitulates experimental results from Figure 5d


allcAMP = logspace(-11, -5);
N = length(allcAMP);
A0B0s = zeros(1,N);
A1B0s = zeros(1,N);
A0B1s = zeros(1,N);
A1B1s = zeros(1,N);
cAMPs = zeros(1,N);

for jj = 1:N

    tspan = [0,1000]; 

    A0B0 = 1e-12;  % to test equilibrium, assume [PKA] << [cAMP]
    A1B0 = 0;
    A0B1 = 0;
    A1B1 = 0;
    cAMP = allcAMP(jj);

    y0 = [A0B0, A1B0, A0B1, A1B1, cAMP];

    options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
    [t,y] = ode23s(@PKAactivationOdes, tspan, y0, options);

    A0B0 = y(:,1);
    A1B0 = y(:,2);
    A0B1 = y(:,3);
    A1B1 = y(:,4);
    cAMP = y(:,5);
    
    A0B0s(jj) = y(end,1);
    A1B0s(jj) = y(end,2);
    A0B1s(jj) = y(end,3);
    A1B1s(jj) = y(end,4);
    cAMPs(jj) = y(end,5);

%     AB_tot = A0B0 + A1B0 + A0B1 + A1B1;    % should be 1.26*10^-4
%     m_AB = mean(AB_tot)
%     v_AB = var(AB_tot)
%     cAMP_tot = A1B0 + A0B1 + 2*A1B1 + cAMP;
%     m_cAMP = mean(cAMP_tot)
%     v_cAMP = var(cAMP_tot)
% 
%     name = {'A0B0', 'A1B0', 'A0B1', 'A1B1', 'cAMP'};
%     
%     figure;
%     for ii = 1:5
%         subplot(2,3,ii);
%         plot(t, y(:,ii), 'LineWidth', 2)
%         xlabel('Time (s)','FontSize',14)
%         ylabel(name(ii), 'Fontsize', 14)
%         title(allcAMP(jj));
%     end
end

figure;
AB_tot = A0B0s + A1B0s + A0B1s + A1B1s;
semilogx(allcAMP, A0B0s ./ AB_tot, 'DisplayName', 'A0B0', 'LineWidth', 2);
hold on;
semilogx(allcAMP, A1B1s ./ AB_tot, 'DisplayName', 'A1B1', 'LineWidth', 2);
semilogx(allcAMP, A1B0s ./ AB_tot, 'DisplayName', 'A1B0', 'LineWidth', 2);
semilogx(allcAMP, A0B1s ./ AB_tot, 'DisplayName', 'A0B1', 'LineWidth', 2);

% Comparing to results from entire holoenzyme (dimer vs. tetramer)
% https://science.sciencemag.org/content/335/6069/712
Ka = 40e-9;
h = 0.81;
semilogx(allcAMP, allcAMP .^ h ./ (Ka^h + allcAMP .^ h), ...
    'DisplayName', 'RIIa(90-400):C activation', 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
Ka = 137e-9;
h = 1.5;
semilogx(allcAMP, allcAMP .^ h ./ (Ka^h + allcAMP .^ h), ...
    'DisplayName', 'RIIa_2:C_2 activation', 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 2);
ylabel('Fraction');
xlabel('[cAMP]');
legend('Location', 'west');
