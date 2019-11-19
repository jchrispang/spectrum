function [nominal_spectra, nominal_normSpectra] = calc_nominal_spectra(f)
%% calc_nominal_spectra.m
%
% Calculates spectra using nominal model parameters
%
% Input: f                    : vector of frequencies
%
% Outputs: nominal_spectra    : structure containining the nominal power spectra 
%                               Possible fields are P0, P1, P2, P3, and PBOLD
%         nominal_normSpectra : structure containining the normalized nominal power spectra 
%                               Possible fields are P0, P1, P2, P3, and PBOLD
%
% Example:
% >> f = linspace(0.01, 1, 1000);
% >> [nominal_spectra, nominal_normSpectra] = calc_nominal_spectra(f);
% >> nominal_spectra.PBOLD  % gives out the BOLD power spectrum
%
% Original: James Pang, QIMR Berghofer Medical Research Institute, 2019

%% main code

params = utils.loadParameters_new;
[nominal_spectra, nominal_normSpectra] = utils.calc_spectra(f, params);
