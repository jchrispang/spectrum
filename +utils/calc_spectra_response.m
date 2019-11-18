function [spectra, BOLD, BOLD_center, f, x, t] = calc_spectra_response(params)
% calculate the different components of the spectra and the BOLD response

%% for spectra
f = linspace(0.001, 1, 10);
w = 2*pi*f;

Y1 = (params.k2 - params.k3);
Y2 = (params.k1 + params.k2);
P = -params.C_z*(Y1 - params.V_0*Y2);
Q = params.C_z*(Y1*(params.eta + params.tau^(-1)) - ...
Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)) + ...
    (params.D/params.rho_f)*(Y1 - params.V_0*Y2));
R = params.C_z*(params.D/params.rho_f)*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)));

constant = (pi/(2*params.v_b^2*params.Gamma));

% numerator
spectra.num = w.^4*P^2 + w.^2*(Q^2 + 2*P*R) + R^2; 
% W
spectra.W = constant*(1./w).*(pi/2 - atan((params.k_z^2*params.v_b^2 - w.^2)./(2*params.Gamma*w)));
% L
spectra.L = 1./((-w.^2 + params.kappa^2/4 + params.w_f^2).^2 + w.^2*params.kappa^2);
% D
spectra.D = 1./(w.^2 + (params.eta + params.tau^(-1)).^2);
% final spectra
spectra.total = spectra.num.*spectra.W.*spectra.L.*spectra.D;

%% for BOLD response
d_max = 10e-3;
t_max = 50;

params.Nkx = 2^9;
params.Nw = 2^9;

x = d_max*(2*(0:params.Nkx-1)/params.Nkx - 1);
t = t_max*(2*(0:params.Nw-1)/params.Nw - 1);
center = dsearchn(x', 0);

kxsamp = (1/mean(diff(x)))*2*pi;
wsamp = (1/mean(diff(t)))*2*pi;

[kx, w] = generate_kw_1D(kxsamp, wsamp, params.Nkx, params.Nw);

T = calcTransFuncs_fromPhi_1D(kx, w, params);
neural_freq = 1;

BOLD = real(freq2coord_1D(T.T_Yphi.*neural_freq, kx, w));
BOLD_center = BOLD(center, :);
