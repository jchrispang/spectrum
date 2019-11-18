function [spectra, normSpectra] = calc_spectra_limits(f, flow, fhigh, params)
% calculate limits of different spectra

w = 2*pi*f;
wlow = 2*pi*flow;
whigh = 2*pi*fhigh;

Y1 = (params.k2 - params.k3);
Y2 = (params.k1 + params.k2);
P = -params.C_z*(Y1 - params.V_0*Y2);
Q = params.C_z*(Y1*(params.eta + params.tau^(-1)) - ...
Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)) + ...
    (params.D/params.rho_f)*(Y1 - params.V_0*Y2));
R = params.C_z*(params.D/params.rho_f)*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)));

constant = 1/(8*pi*params.v_b^2*params.Gamma);
% low-frequency limits
% spectra.low.P0 = R^2*ones(size(f));
% spectra.low.P1 = (1/(4*pi*params.k_z^2*params.v_b^4))*ones(size(f));
% spectra.low.P2 = (1/(params.kappa^2/4 + params.w_f^2)^2)*ones(size(f));
% spectra.low.P3 = (1/(params.eta + params.tau^(-1))^2)*ones(size(f));
% spectra.low.PBOLD = spectra.low.P0.*spectra.low.P1.*spectra.low.P2.*spectra.low.P3;
spectra.low.P0 = R^2*ones(size(flow));
spectra.low.P1 = (1/(4*pi*params.k_z^2*params.v_b^4))*ones(size(flow));
spectra.low.P2 = (1/(params.kappa^2/4 + params.w_f^2)^2)*ones(size(flow));
spectra.low.P3 = (1/(params.eta + params.tau^(-1))^2)*ones(size(flow));
spectra.low.PBOLD = spectra.low.P0.*spectra.low.P1.*spectra.low.P2.*spectra.low.P3;
% high-frequency limits
% spectra.high.P0 = P^2*w.^4;
% spectra.high.P1 = constant*pi*(1./w);
% spectra.high.P2 = 1./(w.^4);
% spectra.high.P3 = 1./(w.^2);
% spectra.high.PBOLD = spectra.high.P0.*spectra.high.P1.*spectra.high.P2.*spectra.high.P3;
spectra.high.P0 = P^2*whigh.^4;
spectra.high.P1 = constant*pi*(1./whigh);
spectra.high.P2 = 1./(whigh.^4);
spectra.high.P3 = 1./(whigh.^2);
spectra.high.PBOLD = spectra.high.P0.*spectra.high.P1.*spectra.high.P2.*spectra.high.P3;

% Normalized spectra
if nargout > 1
%     normSpectra.low.P0 = spectra.low.P0/trapz(w, spectra.low.P0);
%     normSpectra.low.P1 = spectra.low.P1/trapz(w, spectra.low.P1);
%     normSpectra.low.P2 = spectra.low.P2/trapz(w, spectra.low.P2);
%     normSpectra.low.P3 = spectra.low.P3/trapz(w, spectra.low.P3);
%     normSpectra.low.PBOLD = spectra.low.PBOLD/trapz(w, spectra.low.PBOLD);
    normSpectra.low.P0 = spectra.low.P0/trapz(wlow, spectra.low.P0);
    normSpectra.low.P1 = spectra.low.P1/trapz(wlow, spectra.low.P1);
    normSpectra.low.P2 = spectra.low.P2/trapz(wlow, spectra.low.P2);
    normSpectra.low.P3 = spectra.low.P3/trapz(wlow, spectra.low.P3);
    normSpectra.low.PBOLD = spectra.low.PBOLD/trapz(wlow, spectra.low.PBOLD);
    
%     normSpectra.high.P0 = spectra.high.P0/trapz(w, spectra.high.P0);
%     normSpectra.high.P1 = spectra.high.P1/trapz(w, spectra.high.P1);
%     normSpectra.high.P2 = spectra.high.P2/trapz(w, spectra.high.P2);
%     normSpectra.high.P3 = spectra.high.P3/trapz(w, spectra.high.P3);
%     normSpectra.high.PBOLD = spectra.high.PBOLD/trapz(w, spectra.high.PBOLD);
    normSpectra.high.P0 = spectra.high.P0/trapz(whigh, spectra.high.P0);
    normSpectra.high.P1 = spectra.high.P1/trapz(whigh, spectra.high.P1);
    normSpectra.high.P2 = spectra.high.P2/trapz(whigh, spectra.high.P2);
    normSpectra.high.P3 = spectra.high.P3/trapz(whigh, spectra.high.P3);
    normSpectra.high.PBOLD = spectra.high.PBOLD/trapz(whigh, spectra.high.PBOLD);
end