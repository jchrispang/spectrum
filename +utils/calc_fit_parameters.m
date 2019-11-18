function [] = calc_fit_parameters()

data_types = {'He', 'Chen_Boynton', 'Chen_Bandettini'};

params = utils.loadParameters_new;
limits = utils.get_params_limits_new();

opts = optimset('TolFun',1e-6, 'Display', 'off');

init_parameters.(data_types{1}) = [1.2784, 0.34, 0.165];
init_parameters.(data_types{2}) = [1.2784, 0.44, 0.22];
init_parameters.(data_types{3}) = [1.415, 0.8, 0.39];
% fit_parameters = zeros(1,3);
% fit_parameters(1) = params.tau;         % tau
% fit_parameters(2) = params.kappa;       % kappa
% fit_parameters(3) = params.w_f;         % w_f

for i = 1:length(data_types)
    fit_parameters(1) = init_parameters.(data_types{i})(1);
    fit_parameters(2) = init_parameters.(data_types{i})(2);
    fit_parameters(3) = init_parameters.(data_types{i})(3);
    [fit, sse, ~, ~, ~, ~, jacobian] = lsqnonlin('BOLDpower_LSE', ...
                fit_parameters, [limits.tau(1), limits.kappa(1), limits.w_f(1)], ...
                [limits.tau(2), limits.kappa(2), limits.w_f(2)], opts, ...
                data_types{i});
            
    fit_values.(data_types{i}) = fit;
    sse_values.(data_types{i}) = sse;
    jacobian_values.(data_types{i}) = jacobian;
    
    difference = BOLDpower_LSE(fit, data_types{i});
    rmse_values.(data_types{i}) = sqrt(sum(difference.^2));
end

fit_values

save('data/fitted_parameters.mat', 'fit_values', 'sse_values', ...
                                   'jacobian_values', 'rmse_values')
