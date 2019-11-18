function combining_literature_data

% He et al. (2010) data
regions = 21;
He_arrangement = [21, 12, 13, 20, 10, 4, 1, 15, 11, 9, 5, 7, 8, 2, 16, ...
                   3, 18, 17, 19, 14, 6];
               
for i=1:regions
    data = importdata(sprintf('./data/data_He_%d_%d.dat', i, He_arrangement(i)));
    data_He.x(i, :) = data(:, 1)';
    data_He.y(i, :) = data(:, 2)';
end

% 6 attention, 5 default mode, 3 motor, 2 saliency, 2 visual, 3 non-neocortical
data_He.region_names = {'vIPS', 'R TPJ', 'R DLPFC', 'pIPS', 'MT', 'FEF', ...
                        'AG', 'SFG', 'PCC', 'MPF', 'Frontopolar', ...
                        'L S2', 'L motor', 'Broca', ...
                        'R FI', 'dACC', ...
                        'vRetino', 'dRetino', ...
                        'tha', 'R cerebellum', 'HF'};
% data_He.region_names = {'AG', 'Broca', 'dACC', 'FEF', 'Frontopolar', 'HF', 'L S2', ...
%                         'L motor', 'MPF', 'MT', 'PCC', 'R TPJ', 'R DLPFC', 'R cerebellum', ...
%                         'SFG', 'R FI', 'dRetino', 'vRetino', 'tha', 'pIPS', 'vIPS'};
data_He.network = {'attention', 'attention', 'attention', 'attention', 'attention', 'attention', ...
                   'default-mode', 'default-mode', 'default-mode', 'default-mode', 'default-mode', ...
                   'motor', 'motor', 'motor', ...
                   'saliency', 'saliency', ...
                   'visual', 'visual', ...
                   'non-neocortical', 'non-neocortical', 'non-neocortical'};
data_He.mean_x = mean(data_He.x);
data_He.mean_y_all = mean(data_He.y);
data_He.std_y_all = std(data_He.y);
data_He.mean_y_cortical = mean(data_He.y(1:18, :));
data_He.std_y_cortical = std(data_He.y(1:18, :));

% Bandettini (1999) data
data = importdata('./data/data_Bandettini.dat');

data_Bandettini.x = data(:, 1)';
data_Bandettini.y = data(:, 2)';

% Boynton et al. (1996) data
data = importdata('./data/data_Boynton.dat');

data_Boynton.x = data(:, 1)';
data_Boynton.y = data(:, 2)';

% Chen and Tyler (2008) fit of Boynton et al. (1996) data
data = importdata('./data/data_Chen_Boynton.dat');

data_Chen_Boynton.x = data(:, 1)';
data_Chen_Boynton.y = data(:, 2)';
data_Chen_Boynton.interp_x = linspace(min(data_Chen_Boynton.x), max(data_Chen_Boynton.x), 100);
data_Chen_Boynton.interp_y = interp1(data_Chen_Boynton.x, data_Chen_Boynton.y, data_Chen_Boynton.interp_x);

% Chen and Tyler (2008) fit of Bandettini (1999) data
data = importdata('./data/data_Chen_Bandettini.dat');

data_Chen_Bandettini.x = data(:, 1)';
data_Chen_Bandettini.y = data(:, 2)';
data_Chen_Bandettini.interp_x = linspace(min(data_Chen_Bandettini.x), max(data_Chen_Bandettini.x), 100);
data_Chen_Bandettini.interp_y = interp1(data_Chen_Bandettini.x, data_Chen_Bandettini.y, data_Chen_Bandettini.interp_x);


save('./data/data.mat', 'data_He', 'data_Bandettini', 'data_Boynton', ...
                        'data_Chen_Boynton', 'data_Chen_Bandettini')