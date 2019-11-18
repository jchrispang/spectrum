% Function to calculate difference of model and experimental BOLD power spectrum data
% fit_parameters = {tau, kappa, w_f}

% Created by: James Pang 
% Date of creation: Sep 14, 2018

function difference = BOLDpower_LSE(fit_parameters, params, data_type)

%% Load experimental data: either He, Boynton, Bandettini, Chen_Boynton, Chen_Bandettini

if strcmpi(data_type, 'He')
    load('data/data.mat', 'data_He')
    experiment.f = data_He.mean_x;
    experiment.P = data_He.mean_y_cortical;
elseif strcmpi(data_type, 'Boynton')
    load('data/data.mat', 'data_Boynton')
    experiment.f = data_Boynton.x;
    experiment.P = data_Boynton.y;
elseif strcmpi(data_type, 'Bandettini')
    load('data/data.mat', 'data_Bandettini')
    experiment.f = data_Bandettini.x;
    experiment.P = data_Bandettini.y;
elseif strcmpi(data_type, 'Chen_Boynton')
    load('data/data.mat', 'data_Chen_Boynton')
    experiment.f = data_Chen_Boynton.interp_x;
    experiment.P = data_Chen_Boynton.interp_y;
elseif strcmpi(data_type, 'Chen_Bandettini')
    load('data/data.mat', 'data_Chen_Bandettini')
    experiment.f = data_Chen_Bandettini.interp_x;
    experiment.P = data_Chen_Bandettini.interp_y;
end

% normalize power spectrum to maximum value
experiment.P = experiment.P/max(experiment.P);

%% Load model parameter values

params.tau = fit_parameters(1);
params.kappa = fit_parameters(2);
params.w_f = fit_parameters(3);

%% Calculate the model response using the fitting parameters 

[~, normSpectra] = utils.calc_spectra(experiment.f, params);

model.P = normSpectra.BOLD;

% normalize power spectrum to maximum value
model.P = model.P/max(model.P);

% %% Interpolate predicted response to experimental time
% 
% centerResponse_prediction = interp1(PredictedResponse.time, ...
%     real(PredictedResponse.Y_xt(size(PredictedResponse.Y_xt,1)/2+1,:)), ...
%     time_experiment);


%% Calculate difference of model and experiment in logarithmic scale

difference = log10(model.P) - log10(experiment.P);

