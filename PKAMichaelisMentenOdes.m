% Modeling of PKA activation by cAMP binding
% 
% http://www.perkinelmer.co.jp/Portals/0/resource/products_ls/reader/pdf/LC3000-AP-207.pdf
% Figure 4: Michaelis-Menten constants (assume [ATP] >> K_M_ATP = 4.22 µM)

function dydt = PKAMichaelisMentenOdes(~,y,E0)

% System of ODEs describing phosphorylation by PKA
% Inputs:
% y: a 1x2 [initial concentrations]
%           
% Output (output is dydt):
% dydt: a 2x1

S = y(1);
P = y(2);

% rate constants - assumes all cAMP on rates are 1e7 1/M/s
k_cat = 2.7;   % 1/s
K_M = 2.79e-6;  % M

dydt = zeros(2,1);

% dS/dt
dydt(1,:) = -k_cat*E0*S/(S + K_M);
% dP/dt
dydt(2,:) = k_cat*E0*S/(S + K_M);

end