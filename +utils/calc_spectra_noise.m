function out = calc_spectra_noise(params, noise_type)
%% calc_spectra_noise.m
%
% Calculates noise and BOLD power spectra for different noise types 
%
% Inputs: params    : instance of the class utils.loadParameters_new
%        noise_type : type of noise (string)
%                     Possible strings are white, pink, brown, blue
%
% Output: out       : structure of output results 
%                     Possible fields are w, spectra, normSpectra,
%                     noise_spectra, noise_normSpectra
%
% Example:
% >> params = utils.loadParameters_new;
% >> noise_type = 'white';
% >> out = calc_spectra_noise(params, noise_type);
% >> out.spectra  % gives out the BOLD power spectrum
%
% Original: James Pang, QIMR Berghofer Medical Research Institute, 2019

%% main code

% Defining noise
switch noise_type
    case 'white'
        exponent = 0;
    case 'pink'
        exponent = -1;
    case 'brown'
        exponent = -2;
    case 'blue'
        exponent = 1;
end

noise = @(w) w.^exponent;

% Setting up parameters for transfer function calculation
d_max = 8e-3; t_max = 120;
params.Nkx = 2^8; params.Nw = 2^12;

x = d_max*(2*(0:params.Nkx-1)/params.Nkx - 1);
t = t_max*(2*(0:params.Nw-1)/params.Nw - 1);
kxsamp = (1/mean(diff(x)))*2*pi;
wsamp = (1/mean(diff(t)))*2*pi;

[kx, w] = generate_kw_1D(kxsamp, wsamp, params.Nkx, params.Nw);

% Defining positive kx and w
kx_center = dsearchn(kx', 0);
w_center = dsearchn(w', 0);

kx = kx(kx_center:end);
w = w(w_center+1:end);
[~, kxmat] = meshgrid(w, kx);

% Calculating BOLD transfer function
T = calcTransFuncs_fromPhi_1D(kx, w, params);
TBOLD = T.T_Yphi;

% Calculating power from transfer function
Pw_TBOLD = trapz(kx, kxmat.*abs(TBOLD).^2/(2*pi));

% Calculating power from noise
Pw_noise = abs(noise(w)).^1;

% Calculating final spectra
out.w = w;
out.spectra = Pw_TBOLD.*Pw_noise;
out.normSpectra = out.spectra/trapz(w, out.spectra);
out.noise_spectra = noise(w);
out.noise_normSpectra = out.noise_spectra/trapz(w, out.noise_spectra);
