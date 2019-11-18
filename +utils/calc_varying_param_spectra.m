function [varying_param_spectra, varying_param_normSpectra] = calc_varying_param_spectra(f, param_str, number_points)

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
