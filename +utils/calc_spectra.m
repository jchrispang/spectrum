function [spectra, normSpectra] = calc_spectra(f, params)
%% calc_spectra.m
%
% Calculates different components of BOLD power spectra 
%
% Inputs: f           : vector of frequencies
%        params       : instance of the class utils.loadParameters_new
%
% Outputs: spectra    : structure containining the power spectra 
%                       Possible fields are P0, P1, P2, P3, and PBOLD
%         normSpectra : structure containining the normalized nominal power spectra (optional)
%                       Possible fields are P0, P1, P2, P3, and PBOLD
%
% Example:
% >> params = utils.loadParameters_new;
% >> f = linspace(0.01, 1, 1000);
% >> [spectra, normSpectra] = calc_spectra(f, params);
% >> spectra.PBOLD  % gives out the BOLD power spectrum
%
% Original: James Pang, QIMR Berghofer Medical Research Institute, 2019

%% main code

w = 2*pi*f;

Y1 = (params.k2 - params.k3);
Y2 = (params.k1 + params.k2);
P = -params.C_z*(Y1 - params.V_0*Y2);
Q = params.C_z*(Y1*(params.eta + params.tau^(-1)) - ...
Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)) + ...
    (params.D/params.rho_f)*(Y1 - params.V_0*Y2));
R = params.C_z*(params.D/params.rho_f)*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)));

constant = 1/(8*pi*params.v_b^2*params.Gamma);

spectra.P0 = P^2*w.^4 + (Q^2 + 2*P*R)*w.^2 + R^2;
spectra.P1 = constant*(1./w).*(pi/2 - atan((params.k_z^2*params.v_b^2 - w.^2)./(2*params.Gamma*w)));
spectra.P2 = 1./((-w.^2 + params.kappa^2/4 + params.w_f^2).^2 + params.kappa^2*w.^2);
spectra.P3 = 1./(w.^2 + (params.eta + params.tau^(-1)).^2);
spectra.PBOLD = spectra.P0.*spectra.P1.*spectra.P2.*spectra.P3;

% Normalized spectra
if nargout > 1
    normSpectra.P0 = spectra.P0/trapz(w, spectra.P0);
    normSpectra.P1 = spectra.P1/trapz(w, spectra.P1);
    normSpectra.P2 = spectra.P2/trapz(w, spectra.P2);
    normSpectra.P3 = spectra.P3/trapz(w, spectra.P3);
    normSpectra.PBOLD = spectra.PBOLD/trapz(w, spectra.PBOLD);
end
