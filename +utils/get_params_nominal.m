function [nominal, nominal_str] = get_params_nominal()
%% get_params_nominal.m     
%
% Gets nominal values of some parameters of the hemodynamic model based
% from Pang et al. (J. Neurosci. Meth., 2018)
%
% James Pang, University of Sydney, 2018

%%
% limits.D = 591.1;         % effective blood viscosity [kg m^(-3) s^(-1)]

% limits.rho_f = 1062;      % blood mass density [kg m^(-3)]

nominal.tau = 1;            % hemodynamic transit time [s]
                           
% limits.alpha = 0.31;      % Grubb's exponent [unitless] 
                             
nominal.beta = 3.2;         % mean elasticity exponent of cortical vessels [unitless]
             
nominal.kappa = 0.57;       % blood flow signal decay rate [s^(-1)]
                             
nominal.w_f = 0.49;         % natural frequency of flow response [s^(-1)]
                                                   
nominal.L = 3e-3;           % average cortical thickness [m]

nominal.v_b = 2e-3;         % wave sped [m/s]
                      
nominal.Gamma = 0.8;        % wave damping rate [s^(-1)]
                           
% limits.tau_d = 1.2;       % astrocytic delay [s]
                           
% limits.Znorm = [0, 0.9999]; % limits for acos()
                           
if nargout > 1
    nominal_str = {'tau', 'beta', 'kappa', 'w_f', 'L', 'v_b', 'Gamma'};
end