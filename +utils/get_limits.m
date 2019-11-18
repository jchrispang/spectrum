function [limits, limits_str] = get_limits()
%% get_limits.m     
%
% Gets physiological limits of some parameters of the hemodynamic model.
%
% James Pang, University of Sydney, 2017

%%
limits.v_b = [1, 12]*1e-3; % wave sped [m/s]
                           % Aquino et al. (PLoS CB 2012), Aquino et al. (JTB 2014)
                      
limits.Gamma = [0.1, 1];   % wave damping rate [s^(-1)]
                           % Aquino et al. (PLoS CB 2012), Aquino et al. (JTB 2014)
                           
% limits.rho_f = [1054, 1065]; % blood mass density [kg m^(-3)]
%                              % Trudnowski and Rico (Clin Chem, 1974)
                             
limits.alpha = [0.28, 0.56]; % Grubb's exponent [unitless]
                             % Boas and Payne (Phys. Meas., 2009)
                             % From Fung (Biomechanics..., 1993), beta =
                             % [1,5], which translates to alpha = [0.2, 1]
                             
limits.tau = [0,5];%[1, 4];       % hemodynamic transit time [s]
                           % Buxton et al. (NeuroImage, 2004)
                           
limits.L = [1, 4.5]*1e-3;  % average cortical thickness [m]
                           % Fischl and Dale (PNAS, 2000)
                           % more constrained postmortem range [1.8, 3.2] mm
                           
limits.kappa = [0,2];%[0,1];%[0.44, 0.70]; % blood flow signal decay rate [s^(-1)]
                             % Pang et al. (NeuroImage, 2017)
                             
limits.w_f = [0,2];%[0,1];%[0.37, 0.61]; % natural frequency of flow response [s^(-1)]
                           % Pang et al. (NeuroImage, 2017)
                           
limits.Znorm = [0, 0.9999]; % limits for acos()
                           
if nargout > 1
    limits_str = {'v_b', 'Gamma', 'rho_f', 'alpha', 'tau', 'L', 'kappa', 'w_f'};
end
                           