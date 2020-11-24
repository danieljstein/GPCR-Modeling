% Modeling of PKA activation by cAMP binding
% 
% https://www.nature.com/articles/s41467-019-11930-2
% Figure 5: Binding curves with microscopic binding affinities

function dydt = PKAactivationOdes(~,y)

% System of ODEs describing the binding of cAMP to PKA
% Inputs:
% y: a 1x5 [initial concentrations]
%           
% Output (output is dydt):
% dydt: a 5x1

A0B0 = y(1);
A1B0 = y(2);
A0B1 = y(3);
A1B1 = y(4);
cAMP = y(5);

% binding affinities
Kd_A1B0_A0B0 = 17e-9;               % M
Kd_A0B1_A0B0 = 10e-9;               % M
Kd_A1B1_A0B1 = Kd_A1B0_A0B0 / 3;    % M
Kd_A1B1_A1B0 = Kd_A0B1_A0B0 / 3;    % M

% rate constants - assumes all cAMP on rates are 1e7 1/M/s
k_on = 1e6;                 % 1/M/s
k1 = k_on;                  % 1/M/s
k2 = Kd_A1B0_A0B0 * k1;     % 1/s
k3 = k_on;                  % 1/M/s
k4 = Kd_A0B1_A0B0 * k3;     % 1/s
k5 = k_on;                  % 1/M/s
k6 = Kd_A1B1_A1B0 * k5;     % 1/s
k7 = k_on;                  % 1/M/s
k8 = Kd_A1B1_A0B1 * k7;     % 1/s

dydt = zeros(5,1);

% dA0B0/dt
dydt(1,:) = -k1*A0B0*cAMP + k2*A1B0 - k3*A0B0*cAMP + k4*A0B1;
% dA1B0/dt
dydt(2,:) = k1*A0B0*cAMP - k2*A1B0 - k5*A1B0*cAMP + k6*A1B1;
% dA0B1/dt
dydt(3,:) = k3*A0B0*cAMP - k4*A0B1 - k7*A0B1*cAMP + k8*A1B1;
% dA1B1/dt
dydt(4,:) = k5*A1B0*cAMP - k6*A1B1 + k7*A0B1*cAMP - k8*A1B1;
% dcAMP/dt
dydt(5,:) = - k1*A0B0*cAMP + k2*A1B0 - k3*A0B0*cAMP + k4*A0B1 ...
            - k5*A1B0*cAMP + k6*A1B1 - k7*A0B1*cAMP + k8*A1B1;

end 