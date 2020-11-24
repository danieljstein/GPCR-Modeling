% Kinetic modeling of adenylyl cyclase
% 
% https://www.jbc.org/content/272/44/27787.full.pdf
% Figure 5: Model for catalytic mechanism
% Table II: Rate constants to fit experimental data

function dydt = ACodes(~,y)

% System of ODEs describing the AC kinetics
% Inputs:
% y: a 1x8 [initial concentrations]
%           
% Output (output is dydt):
% dydt: a 8x1

E = y(1);   % E
S = y(2);   % ATP
ES = y(3);  % E-ATP
EQP = y(4); % E-cAMP-PP
EQ = y(5);  % E-cAMP
EP = y(6);  % E-PP
Q = y(7);   % cAMP
P = y(8);   % PP

% rate constants
k1 = 2.62e5;    % M/s (should this be 1/M/s?)
k2 = 89.5;      % 1/s
k3 = 59;        % 1/s
k4 = 2.6;       % 1/s
k5 = 0.8;       % 1/s
k6 = 2.78e3;    % M/s
k7 = 1060;      % 1/s
k8 = 1.11e5;    % M/s
k9 = 0.39;      % 1/s
k10 = 142;      % M/s
k11 = 56;       % 1/s
k12 = 3.54e5;   % M/s


dydt = zeros(8,1);

% dE/dt
dydt(1,:) = -k1*E*S + k2*ES + k9*EQ - k10*E*Q + k11*EP - k12*E*P;
% dS/dt
dydt(2,:) = -k1*E*S + k2*ES;
% dES/dt
dydt(3,:) = k1*E*S - k2*ES - k3*ES + k4*EQP;
% dEQP/dt
dydt(4,:) = k3*ES - k4*EQP - k5*EQP + k6*EQ*P - k7*EQP + k8*EP*Q;
% dEQ/dt
dydt(5,:) = k5*EQP - k6*EQ*P - k9*EQ + k10*E*Q;
% dEP/dt
dydt(6,:) = k7*EQP - k8*EP*Q - k11*EP + k12*E*P;
% dQ/dt
dydt(7,:) = k7*EQP - k8*EP*Q + k9*EQ - k10*E*Q;
% dP/dt
dydt(8,:) = k5*EQP - k6*EQ*P + k11*EP - k12*E*P;

end 