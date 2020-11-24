function dydt = PolymersomeOdes(~,y)

% System of ODEs describing the polymersome signaling
% Incorporates GPCR signaling, AC activity, and PKA activation
% Inputs:
% y: a 1x12 [initial concentrations]

% Output (output is dydt):
% dydt: a 12x1

yc = num2cell(y);

[L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp] = yc{1:17};
y_GPCR = [L, R, G, Rp, AR, ARp, RG, RpG, ARG, ARpG, Ga_GTP, Gby, Ga_GDP, GTP, GDP, AC, Ga_GTP_ACp];
dydt_GPCR = GPCRodes4(nan, y_GPCR);

[E, S, ES, EQP, EQ, EP, Q, P] = yc{18:25}; % should have E + ES + EQP + EP = Ga_GTP_ACp
E_tot = E + ES + EQP + EQ + EP;
y_AC = [E, S, ES, EQP, EQ, EP, Q, P];
dydt_AC = ACodes(nan, y_AC);

[A0B0, A1B0, A0B1, A1B1, cAMP] = yc{26:30};
y_PKA = [A0B0, A1B0, A0B1, A1B1, cAMP];
dydt_PKA = PKAactivationOdes2(nan, y_PKA);

[S, P] = yc{31:32};
y_PKA_MM = [S, P];
dydt_PKA_MM = PKAMichaelisMentenOdes(nan,y_PKA_MM,A1B1); % assumes that all PKA bound to 2 cAMPs are active

dydt = zeros(32:1);

dydt(1:17,:) = dydt_GPCR;
dydt(18:25,:) = dydt_AC;
dydt(26:30,:) = dydt_PKA;
dydt(31:32,:) = dydt_PKA_MM;

% dydt(18,:) = dydt(18,:) + dydt_GPCR(17,:);

if E_tot > 0
    dAC = dydt_GPCR(17,:)*y([18,20,21,22,23]) / E_tot;
    dydt([18,20,21,22,23],:) = dydt([18,20,21,22,23],:) + dAC;
    dydt(19,:) = dydt(19,:) - dAC(2,:);             % S = ATP
    dydt(24,:) = dydt(24,:) - dAC(3,:) - dAC(4,:);  % Q = cAMP
    dydt(25,:) = dydt(25,:) - dAC(3,:) - dAC(5,:);  % P = PP
else
    dydt(18,:) = dydt(18,:) + dydt_GPCR(17,:);
end
dydt(24,:) = dydt(24,:) + dydt_PKA(5,:);        % cAMP
dydt(30,:) = dydt(24,:);                        % cAMP
dydt(19,:) = dydt(19,:) + dydt_PKA_MM(1,:);     % ATP

end 