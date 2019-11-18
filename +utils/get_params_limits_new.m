function [limits, limits_str] = get_params_limits_new()
%% get_params_limits_new.m     
%
% Gets physiological limits of some parameters of the hemodynamic model.
%
% James Pang, University of Sydney, 2018

%%
% limits.D = [106, 850];     % effective blood viscosity [kg m^(-3) s^(-1)]

% limits.rho_f = [1054, 1065]; % blood mass density [kg m^(-3)]
%                              % Trudnowski and Rico (Clin Chem, 1974)

limits.tau = [1, 4];       % hemodynamic transit time [s]
                           % Buxton et al. (NeuroImage, 2004)
                           
% limits.alpha = [0.28, 0.56]; % Grubb's exponent [unitless]
%                              % Boas and Payne (Phys. Meas., 2009)
%                              % From Fung (Biomechanics..., 1993), beta =
%                              % [1,5], which translates to alpha = [0.2, 1] 
                             
limits.beta = [1/0.56, 1/0.28]; % mean elasticity exponent of cortical vessels [unitless]
             
limits.kappa = [0.1, 1];%[0.44, 0.70]; % blood flow signal decay rate [s^(-1)]
                                       % Pang et al. (NeuroImage, 2017)
                             
limits.w_f = [0.1, 1];%[0.37, 0.61]; % natural frequency of flow response [s^(-1)]
                           % Pang et al. (NeuroImage, 2017)
                                                   
limits.L = [1, 4.5]*1e-3;  % average cortical thickness [m]
                           % Fischl and Dale (PNAS, 2000)
                           % more constrained postmortem range [1.8, 3.2] mm

limits.v_b = [1, 12]*1e-3; % wave sped [m/s]
                           % Aquino et al. (PLoS CB 2012), Aquino et al. (JTB 2014)
                      
limits.Gamma = [0.1, 1];   % wave damping rate [s^(-1)]
                           % Aquino et al. (PLoS CB, 2012), Aquino et al. (JTB 2014)
                           
% limits.tau_d = [0.2, 2.4]; % astrocytic delay [s]
%                            % Pang et al. (NeuroImage, 2017]
                           
% limits.Znorm = [0, 0.9999]; % limits for acos()
                           
if nargout > 1
    limits_str = {'tau', 'beta', 'kappa', 'w_f', 'L', 'v_b', 'Gamma'};
end
                           