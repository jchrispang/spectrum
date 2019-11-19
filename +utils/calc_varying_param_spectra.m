function [varying_param_spectra, varying_param_normSpectra] = calc_varying_param_spectra(f, param_str, number_points)
%% calc_varying_param_spectra.m
%
% Recalculates spectra for number_points values of param_str 
%
% Inputs: f                         : vector of frequencies
%        param_str                  : string of parameter name (see utils.loadParameters_new for possible strings)
%        number_points              : number of incremental values of param_str
%                                     from its minimum to maximum to be used
%
% Outputs: varying_param_spectra    : cell of structure containining the power spectra 
%                                     Possible struct fields are P0, P1, P2, P3, and PBOLD
%         varying_param_normSpectra : cell of structure containining the normalized nominal power spectra (optional)
%                                     Possible struct fields are P0, P1, P2, P3, and PBOLD
%
% Example:
% >> f = linspace(0.01, 1, 1000);
% >> param_str = 'tau';
% >> number_points = 10;
% >> [varying_param_spectra, varying_param_normSpectra] = calc_varying_param_spectra(f, param_str, number_points);
% >> varying_param_spectra{1}.PBOLD  % gives out the BOLD power spectrum using the minimum value of param_str
% >> varying_param_spectra{number_points}.PBOLD  % gives out the BOLD power spectrum using the maximum value of param_str
%
% Original: James Pang, QIMR Berghofer Medical Research Institute, 2019

%% main code

% forces a minimum of two points
if number_points < 2
    number_points = 2;
end

params = utils.loadParameters_new;
limits = utils.get_params_limits_new();

param_values = linspace(limits.(param_str)(1), limits.(param_str)(2), number_points);
for j = 1:number_points
    params.(param_str) = param_values(j);
    [varying_param_spectra{j}, varying_param_normSpectra{j}] = utils.calc_spectra(f, params);
end
