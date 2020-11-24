function dydt = GPCRodes4(~,y,Kd)

% System of ODEs describing the GPCR signaling
% based on Shea et al. 2000 (https://pubmed.ncbi.nlm.nih.gov/11205879/)
% with constants taken from Woodroffe et al. 2009
%
% Also adds binding to adenylyl cyclase, which serves as a GAP
%
% Inputs:
% y: a 1x17 [initial concentrations]

% Output (output is dydt):
% dydt: a 17x1

y = num2cell(y);
[L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp] = y{:};

% parameters:
V_polymersome = 0.5e-15; % in liters
N_av = 6.022e23; % avogadro's number

% rate constants
zeta = 1000;
mu = 2;
nu = 1;
k1 = 1;                 % 1/s   (k_p) -> K_P = 1e-3
k2 = 1000;              % 1/s   (k_q) -> K_P = 1e-3
k3 = 8.4e7;             % 1/M/s (k_f) -> 1/K_A = 4.4e-09 M
k4 = 0.37;              % 1/s   (k_r) -> 1/K_A = 4.4e-09 M
k5 = k3;                % 1/M/s (k_f * zeta+) -> 1/(zeta*K_A) = 4.3e-12 M
k6 = k4 / zeta;         % 1/s   (k_r * zeta-) -> 1/(zeta*K_A) = 4.3e-12 M
k7 = k1;                % 1/s   (k_p * zeta+) -> zeta*K_P = 1
k8 = k2 / zeta;         % 1/s   (k_p * zeta-) -> zeta*K_P = 1
k9 = 3.6e7;             % 1/M/s (k+) -> 1/K_G = 8.3e-11 M
k10 = 3e-3;             % 1/s   (k-) -> 1/K_G = 8.3e-11 M
k11 = k9;               % 1/M/s (k+ * mu+) -> 1/(mu*K_G) = 4.2e-11 M
k12 = k10 / mu;         % 1/s   (k- * mu-) -> 1/(mu*K_G) = 4.2e-11 M
k13 = k1;               % 1/s   (k_p * mu+) -> mu*K_P = 2e-3
k14 = k2 / mu;          % 1/s   (k_q * mu-) -> mu*K_P = 2e-3
k15 = k9;               % 1/M/s (k+ * nu+) -> 1/(nu*K_G) = 8.3e-11 M
k16 = k10 / nu;         % 1/s   (k- * nu-) -> 1/(nu*K_G) = 8.3e-11 M
k17 = k9;               % 1/M/s (k+ * mu+ * nu+) -> 1/(mu*nu*K_G) = 4.2e-11 M
k18 = k10 / mu / nu;    % 1/s   (k- * mu- * nu-) -> 1/(mu*nu*K_G) = 4.2e-11 M
k19 = k1;               % 1/s   (k_p * zeta+ * mu+) -> zeta*mu*K_P = 2
k20 = k2 / zeta / mu;   % 1/s   (k_q * zeta- * mu-) -> zeta*mu*K_P = 2
k21 = k3;               % 1/M/s (k_f * zeta+ * nu+) -> 1/(zeta*nu*K_A) = 4.3e-12 M
k22 = k4 / zeta / nu;   % 1/s   (k_r * zeta- * nu-) -> 1/(zeta*nu*K_A) = 4.3e-12 M
k23 = k3;               % 1/M/s (k3 * nu+) -> 1/(nu*K_A) = 4.4e-9 M
k24 = k4 / nu;          % 1/s   (k4 * nu-) -> 1/(nu*K_A) = 4.4e-9 M
kact = 1;               % 1/s   (kb+, assumed irreversible as in Shea et al.)
                        %       note: depends on GTP:GDP
kGTP_slow = 0.1;        % 1/s   (kg+)
kGTP_GAP = 5;           % 1/s   (from RGS modeling in Shea et al. Fig. 4)
k_assoc = 1.2e10;       % 1/M/s (kb+, assumed irreversible as in Shea et al.)
                        %       assume fast reassociation of Ga_GDP and Gby
kAC_on = 1e6;          % 1/M/s (typical protein-protein interaction)
kAC_off = 0.05;        % 1/s   (typical protein-protein interaction)

dydt = zeros(12,1);

% L
dydt(1,:) = - k3*R*L + k4*AR - k5*Rp*L + k6*ARp - k21*RpG*L + k22*ARpG - k23*RG*L + k24*ARG;
% R
dydt(2,:) = - k1*R + k2*Rp - k3*R*L + k4*AR - k9*R*G + k10*RG;
% G
dydt(3,:) = - k9*R*G + k10*RG - k11*Rp*G + k12*RpG - k15*AR*G + k16*ARG - k17*ARp*G + k18*ARpG + k_assoc*Ga_GDP*Gby;
% Rp
dydt(4,:) = k1*R - k2*Rp - k5*Rp*L + k6*ARp - k11*Rp*G + k12*RpG + kact*RpG;
% AR
dydt(5,:) = k3*R*L - k4*AR - k7*AR + k8*ARp - k15*AR*G + k16*ARG;
% ARp
dydt(6,:) = k5*Rp*L - k6*ARp + k7*AR - k8*ARp - k17*ARp*G + k18*ARpG + kact*ARpG;
% RG
dydt(7,:) = k9*R*G - k10*RG - k13*RG + k14*RpG - k23*RG*L + k24*ARG;
% RpG
dydt(8,:) = k11*Rp*G - k12*RpG + k13*RG - k14*RpG - k21*RpG*L + k22*ARpG - kact*RpG;
% ARG
dydt(9,:) = k15*AR*G - k16*ARG - k19*ARG + k20*ARpG + k23*RG*L - k24*ARG;
% ARpG
dydt(10,:) = k17*ARp*G - k18*ARpG + k19*ARG - k20*ARpG + k21*RpG*L - k22*ARpG - kact*ARpG;
% Ga_GTP
dydt(11,:) = kact*RpG + kact*ARpG - kGTP_slow*Ga_GTP ...
             - kAC_on*AC*Ga_GTP + kAC_off*Ga_GTP_ACp;
% Gby
dydt(12,:) = kact*RpG + kact*ARpG - k_assoc*Ga_GDP*Gby;
% Ga_GDP
dydt(13,:) = kGTP_slow*Ga_GTP - k_assoc*Ga_GDP*Gby + kGTP_GAP*Ga_GTP_ACp;
% GTP
dydt(14,:) = -kact*RpG - kact*ARpG;
% GDP
dydt(15,:) = kact*RpG + kact*ARpG;
% AC
dydt(16,:) = - kAC_on*AC*Ga_GTP + kAC_off*Ga_GTP_ACp + kGTP_GAP*Ga_GTP_ACp;
% Ga_GTP_ACp
dydt(17,:) = kAC_on*AC*Ga_GTP - kAC_off*Ga_GTP_ACp - kGTP_GAP*Ga_GTP_ACp;

end 