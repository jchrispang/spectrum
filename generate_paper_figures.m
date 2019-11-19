function generate_paper_figures
issave = 1;
folder_name = 'figures';

if issave
    if ~exist(folder_name, 'dir')
        mkdir(folder_name)
    end
end

figure1(issave,folder_name,'Figure1')
figure4(issave,folder_name,'Figure4')
figure5(issave,folder_name,'Figure5')
figure6(issave,folder_name,'Figure6')
figure7(issave,folder_name,'Figure7')
figure8(issave,folder_name,'Figure8')
figure9(issave,folder_name,'Figure9')
figure10(issave,folder_name,'Figure10')
close all

function figure1(issave,folder_name,file_name)
% Figure 1: Experimental data

load data/data.mat

% He et al. (2010) data
% covers different networks: 6 attention, 5 default-mode, 3 motor,
% 2 saliency, 2 visual, 3 non-neocortical regions
% regions = 21;
regions_cortical = 18;
x = data_He.mean_x;

fig = figure('Position', [200, 200, 750, 350]);
axes('Position', [0.1 0.15 0.38 0.75])
loglog(NaN, NaN)
hold on;
for i=1:regions_cortical
    y = data_He.y(i,:);
    loglog(x, y, '.-', 'linewidth', 1)
end
hold off;

text(0.06, 8, '$s \approx 0.5$ to $1.2$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')

annotation('textbox', [0.14, 0.18, 0.1, 0.1], 'string', 'He at al. [16]', ...
           'fontsize', 15, 'fitboxtotext', 'on', 'margin', 2, ...
           'horizontalalignment', 'center', 'verticalalignment', 'middle', 'linestyle', 'none')

set(gca, 'Xlim', [0.01, 0.2], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2, 1, 5, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1, 5, 10, 100}, ...
         'Ytick', [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ...
         'Ylim', [0.8, 20], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);
xlabel('$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')

% Boynton et al. (1996), Bandettini (1999), and Chen and Tyler (2008)
axes('Position', [0.58 0.15 0.38 0.75])
loglog(NaN, NaN, 'HandleVisibility', 'off')
hold on;
loglog(data_Boynton.x, data_Boynton.y, 'k^', 'markersize', 12)
loglog(data_Bandettini.x, data_Bandettini.y, 'ko', 'markersize', 12)
loglog(data_Chen_Boynton.x, data_Chen_Boynton.y, 'k--', 'linewidth', 2)
loglog(data_Chen_Bandettini.x, data_Chen_Bandettini.y, 'k--', 'linewidth', 2)
hold off;
leg = legend('Boynton et al. [22]', 'Bandettini [23]', 'Chen and Tyler [24]');

set(gca, 'Xlim', [0.01, 0.2], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2, 1, 5, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1, 5, 10, 100}, ...
         'Ylim', [0.08, 1.1], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);
set(leg, 'fontSize', 15, 'location', 'southwest', 'box', 'off');
xlabel('$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')

text(0.045, 0.45, '$s \approx 1.2$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')
text(0.07, 0.2, '$s \approx 2.9$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')

annotation('textbox', [0.03, 0.88, 0.1, 0.1], 'string', '(a)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.51, 0.88, 0.1, 0.1], 'string', '(b)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')

if issave
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure4(issave,folder_name,file_name)
% Figure 4: Spatiotemporal and temporal profile of BOLD for nominal
% parameters

h.p = utils.loadParameters_new;
h = calc_BOLD(h);
h.BOLD_center = h.BOLD(h.center, :);
clim = [min(h.BOLD), max(h.BOLD)];
h.clims = [-max(clim), max(clim)];

h.fig = figure('Position', [200, 200, 650, 250]);
h.ax = axes('Parent', h.fig, 'Position', [0.11 0.19 0.37 0.75]);
prepare_figBOLDspatiotemporal(h)
hold on;
plot([0 0], [h.t(1) h.t(end)], 'k--', 'linewidth', 1.5)
hold off;
set(h.ax.Colorbar.Label, 'String', 'Spatiotemporal response, $Y(x,t)$', ...
    'fontsize', 15, 'interpreter', 'latex');

h.ax = axes('Parent', h.fig, 'Position', [0.69 0.19 0.29 0.75]);
prepare_figBOLDtemporal(h)
ylabel(h.ax, 'Temporal response, $Y(x=0,t)$', 'fontsize', 15, 'interpreter', 'latex');

annotation('textbox', [0.01, 0.92, 0.1, 0.1], 'string', '(a)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.55, 0.92, 0.1, 0.1], 'string', '(b)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure5(issave,folder_name,file_name)
% Figure 5: Spectra using nominal parameters

h.fmin = 0.01;
h.fmax = 1;
h.f = linspace(h.fmin, h.fmax, 1000);
[~, h.normSpectra] = utils.calc_nominal_spectra(h.f);
linecolors = [0 0 0; 215 48 31; 252 141 89; 253 204 138; 254 240 217]/255;
linetypes = {'-.', '-', '-', '-', '-'};
linewidths = 2.5*ones(1,5);

what_spectras = {'PBOLD', 'P0', 'P1', 'P2', 'P3'};
h.fig = figure;
h.ax = axes('Parent', h.fig);
hsetup = loglog(h.ax, NaN, NaN, 'HandleVisibility', 'off');
hold on;
for j = 1:length(what_spectras)
    normalized_spectra = h.normSpectra.(what_spectras{j});
    ylims(j, :) = [min(normalized_spectra), max(normalized_spectra)];
    loglog(h.ax, h.f, normalized_spectra, linetypes{j}, 'linewidth', linewidths(j), ...
           'color', linecolors(j,:));
end
plot([0.3 1], [min(ylims(:,1))*0.8, min(ylims(:,1))*0.8], 'k-', 'linewidth', 3)
plot([0.01 0.02], [min(ylims(:,1))*0.8, min(ylims(:,1))*0.8], 'k-', 'linewidth', 3)
text(0.012, min(ylims(:,1))*1.2, 'LOW', 'fontsize', 15, 'fontweight', 'b')
text(0.44, min(ylims(:,1))*1.2, 'HIGH', 'fontsize', 15, 'fontweight', 'b')
hold off;
leg = legend('$P_{\rm BOLD}$', '$P_0$', '$P_1$', '$P_2$', '$P_3$');

set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2, 1, 5, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1, 5, 10, 100}, ...
         'Ylim', [min(ylims(:,1))*0.8, max(ylims(:,2))*1.2], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);
set(leg, 'fontSize', 15, 'location', 'west', 'interpreter', 'latex', 'box', 'off');
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, '$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure6(issave,folder_name,file_name)
% Figure 6: Autocorrelation

h.p = utils.loadParameters_new;

d_max = 2e-3;
t_max = 120;

h.p.Nkx = 2^8;
h.p.Nw = 2^12;

h.x = d_max*(2*(0:h.p.Nkx-1)/h.p.Nkx - 1);
h.t = t_max*(2*(0:h.p.Nw-1)/h.p.Nw - 1);
h.center = dsearchn(h.t', 0);

kxsamp = (1/mean(diff(h.x)))*2*pi;
wsamp = (1/mean(diff(h.t)))*2*pi;

[~, w] = generate_kw_1D(kxsamp, wsamp, h.p.Nkx, h.p.Nw);

w(h.center) = eps;
[nominal_spectra, ~] = utils.calc_nominal_spectra(w/(2*pi));

wvals = (-1).^(1:length(w));
wM = wvals;
out = wM.*ifft(wM.*abs(nominal_spectra.PBOLD).^2);

autocorr = abs(out(h.center:end));
norm_autocorr = autocorr/autocorr(1);
t = h.t(h.center:end);

% Fitting an exponential decay 
fo = fitoptions('Method', 'NonlinearLeastSquares');
ft = fittype('a*exp(-x/b)+c', 'options', fo);

[curve, gof] = fit(t', norm_autocorr', ft);

h.fig = figure;
h.ax = axes('Parent', h.fig);
plot(t, norm_autocorr, 'b-', 'linewidth', 1.5)
hold on;
plot(t, curve(t), 'k--', 'linewidth', 1.5)
hold off;
leg = legend(h.ax, 'model', 'fit');
xlabel(h.ax, '$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, 'normalized autocorrelation, $C(t)/C(0)$', 'fontsize', 15, 'interpreter', 'latex')
set(leg, 'fontSize', 15, 'location', 'northeast', 'box', 'off');
set(h.ax, 'Xtick', 0:20:120, ...
    'fontsize', 15, 'ticklength', [0.02, 0.02]); 
annotation('textarrow',[0.5 0.2],[0.5 0.7],'String','$\frac{C(t)}{C(t=0)}\sim{\rm e}^{-\frac{t}{\tau_{s}}}$', ...
    'interpreter', 'latex', 'fontsize', 30, 'color', 'k')

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure7(issave,folder_name,file_name)
% Figure 7: Spectra using different tau, kappa, and w_f

h.fmin = 0.01;
h.fmax = 0.2;
h.f = linspace(h.fmin, h.fmax, 1000);
h.number_points = 3;
h.linecolors = [0 0 0; 0 0 1; 1 0 0];
h.linetypes = {'-', '--', '-.'};
h.linewidths = 1.5*ones(1,3);
h.titles = {'$P_{\rm BOLD}(f)$', '$P_0(f)$', '$P_1(f)$', '$P_2(f)$', '$P_3(f)$'};
params_str = {'tau', 'kappa', 'w_f'};
params_names = {'$\tau$', '$\kappa$', '$\omega_f$'};

init_x = 0.12; init_y = 0.68; spacing_x = 1/5.5; spacing_y = 1/5.3;
width = 1/6.8; height = 1/6.2;

h.fig = figure('Position', [200, 200, 800, 500]);
for i = 1:length(params_str)
    [~, varying_param_normSpectra] = utils.calc_varying_param_spectra(h.f, params_str{i}, h.number_points);
    for j = 1:5
        if j==1
            what_spectra = 'PBOLD';
        else
            what_spectra = ['P', num2str(j-2)];
        end
        
        h.ax = axes('Parent', h.fig, ...
                    'Position', [init_x+(j-1)*spacing_x+(j==1)*(-0.04) init_y-(i-1)*spacing_y width height]);
        prepare_figVaryingSpectra(h, varying_param_normSpectra, what_spectra)
        set(h.ax, 'Yticklabel', {}, 'ticklength', [0.03, 0.03])
        
        if i==1
            title(h.ax, h.titles{j}, 'fontsize', 18, 'fontweight', 'b', ...
                                     'interpreter', 'latex')
        end
        if j==1
            annotation('textbox', [0, init_y-(i-1)*spacing_y, 0, height], ...
                       'string', params_names{i}, 'edgecolor', 'none', ...
                       'verticalalignment', 'middle', 'fontsize', 24, ...
                       'fontweight', 'b', 'interpreter', 'latex')
        end
        if i~=length(params_str)
            set(h.ax, 'Xticklabel', {});
        else
            xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
        end
        if i==1 && j==1
            leg = legend('minimum', 'middle', 'maximum');
            set(leg, 'fontSize', 20, 'position', [0.52 0.96 0.01 0.01], ...
                'box', 'off', 'orientation', 'horizontal');
        end
        if prod(varying_param_normSpectra{1}.(what_spectra)./varying_param_normSpectra{2}.(what_spectra))==1
            chH = get(gca,'Children');
            set(gca,'Children',[chH(2); chH(1); chH(end)])
        end
        if i==length(params_str) && j==1
            annotation('rectangle', ...
                       'position', [init_x+(j-1)*spacing_x+(j==1)*(-0.04)-0.03 init_y-(i-1)*spacing_y-0.1 width*1.41 height*4.4], ...
                       'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 2);
        end
    end
end


if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure8(issave,folder_name,file_name)
% Figure 8: Spectra using fitted parameters 

load data/data.mat

h.fmin = 0.003;
h.fmax = 0.2;
h.f = linspace(h.fmin, h.fmax, 1000);

init_x = 0.08; init_y = 0.15; spacing_x = 0.32; 
width = 0.26; height = 0.75;

h.fig = figure('Position', [200, 200, 1000, 350]);

% tuned parameters for He et al. (2010) data
params = utils.loadParameters_new;
params.tau = 1.1;
params.kappa = 0.8;
params.w_f = 0.19;

[~, normSpectra] = utils.calc_spectra(h.f, params);
h.linewidth = 2;
h.linetype = 'k-';

h.ax = axes('Parent', h.fig, ...
             'Position', [init_x+(1-1)*spacing_x init_y width height]);
prepare_figspectrum(h, normSpectra.PBOLD/max(normSpectra.PBOLD));
hold on;
errorbar(h.ax, data_He.mean_x, data_He.mean_y_cortical/max(data_He.mean_y_cortical), ...
         data_He.std_y_cortical/max(data_He.mean_y_cortical), 'bs', ...
         'markersize', 8);
hold off;
set(h.ax, 'ylim', [0.06, 1.3]);
leg = legend(h.ax, 'model', 'He et al. [16]');

xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, '$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')
set(leg, 'fontSize', 15, 'location', 'southwest', 'box', 'off');
text(0.055, 0.16, '$s \approx 3$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')
text(0.027, 0.35, '$s \approx 0.9$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'b', 'interpreter', 'latex')
% title('averaged data', 'fontsize', 15) 

% tuned parameters for Boynton et al. (1996) data
params = utils.loadParameters_new;
params.tau = 1.8;
params.kappa = 1;
params.w_f = 0.36;

[~, normSpectra] = utils.calc_spectra(h.f, params);
h.linewidth = 2;
h.linetype = 'k-';

h.ax = axes('Parent', h.fig, ...
             'Position', [init_x+(2-1)*spacing_x init_y width height]);
prepare_figspectrum(h, normSpectra.PBOLD/max(normSpectra.PBOLD));
hold on;
loglog(h.ax, data_Boynton.x, data_Boynton.y/max(data_Boynton.y), 'b^', 'markersize', 8)
loglog(h.ax, data_Chen_Boynton.x, data_Chen_Boynton.y/max(data_Chen_Boynton.y), 'b--', 'linewidth', 2)
hold off;
leg = legend(h.ax, 'model', 'Boynton et al. [22]', 'Chen and Tyler [24]');

set(h.ax, 'ylim', [0.06, 1.3]);
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
set(leg, 'fontSize', 15, 'location', 'southwest', 'box', 'off');
text(0.055, 0.2, '$s \approx 3$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')
text(0.07, 0.75, '$s \approx 1.2$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'b', 'interpreter', 'latex')
% title('averaged data', 'fontsize', 15) 

% tuned parameters for Bandettini (1999) data
params = utils.loadParameters_new;
params.tau = 1.415;
params.kappa = 0.8;
params.w_f = 0.39;

[~, normSpectra] = utils.calc_spectra(h.f, params);
h.linewidth = 2;
h.linetype = 'k-';

h.ax = axes('Parent', h.fig, ...
             'Position', [init_x+(3-1)*spacing_x init_y width height]);
prepare_figspectrum(h, normSpectra.PBOLD/max(normSpectra.PBOLD));
hold on;
loglog(h.ax, data_Bandettini.x, data_Bandettini.y/max(data_Bandettini.y), 'bo', 'markersize', 8)
loglog(h.ax, data_Chen_Bandettini.x, data_Chen_Bandettini.y/max(data_Chen_Bandettini.y), 'b--', 'linewidth', 2)
hold off;
leg = legend(h.ax, 'model', 'Bandettini [23]', 'Chen and Tyler [24]');

set(h.ax, 'ylim', [0.06, 1.3]);
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
set(leg, 'fontSize', 15, 'location', 'southwest', 'box', 'off');
text(0.065, 0.1819, '$s \approx 3$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')
text(0.045, 0.25, '$s \approx 2.9$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'b', 'interpreter', 'latex')
% title('nonaveraged data', 'fontsize', 15) 

annotation('textbox', [init_x+(1-1)*spacing_x-0.04 init_y+height-0.02, 0.1, 0.1], 'string', '(a)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [init_x+(2-1)*spacing_x-0.04 init_y+height-0.02, 0.1, 0.1], 'string', '(b)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [init_x+(3-1)*spacing_x-0.04 init_y+height-0.02, 0.1, 0.1], 'string', '(c)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure9(issave,folder_name,file_name)
% Figure 9: Creating artificial power laws
% Averaging of multiple voxels, subjects and hemisphere

h.fig = figure;

%%%%%%% Artificial 1: Averaging from random parameters
h.fmin = 0.01;
h.fmax = 0.2;
h.f = linspace(h.fmin, h.fmax, 1000);
params = utils.loadParameters_new;
limits = utils.get_params_limits_new();
what_spectra = 'PBOLD';
param_str = 'w_f'; 

h.number_points = 40;
rng(500000)
param_values = limits.(param_str)(1) + ...
               (params.(param_str) - limits.(param_str)(1))*rand(1, h.number_points);

h.linecolors = cbrewer('seq', 'YlGnBu', h.number_points+1, 'pchip');
h.linetypes = repmat({'-'}, 1, h.number_points+1);
h.linewidths = ones(1, h.number_points+1);
all_spectra = zeros(h.number_points, length(h.f));
all_normSpectra = zeros(h.number_points, length(h.f));
for i=1:h.number_points
    params.(param_str) = param_values(i);
    params.kappa = 2*params.w_f + 0.01;
    [temp_spectra, temp_normSpectra] = utils.calc_spectra(h.f, params);
    all_spectra(i,:) = temp_spectra.(what_spectra);
    all_normSpectra(i,:) = temp_normSpectra.(what_spectra);
end

h.linewidth = 2;
h.linetype = 'k-';

% main figure  
h.ax = axes('Parent', h.fig);
mean_spectra = mean(all_spectra)/trapz(2*pi*h.f, mean(all_spectra));
prepare_figspectrum(h, mean_spectra);
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, '$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')
text(0.09, 0.53, '$s \approx 1.8$', 'fontsize', 15, 'fontweight', 'b', ...
    'color', 'k', 'interpreter', 'latex')

f = h.f;

% inset
h.ax = axes('Parent', h.fig, 'Position', [0.2 0.2 0.35 0.35]);
loglog(h.ax, NaN, NaN, 'HandleVisibility','off')
hold on
for j=1:h.number_points
    normalized_spectra = all_normSpectra(j,:);
    ylims(j,:) = [min(normalized_spectra), max(normalized_spectra)];
    loglog(h.ax, h.f, normalized_spectra, h.linetypes{j+1}, 'linewidth', h.linewidths(j+1), ...
           'color', h.linecolors(j,:));
end
hold off
set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 1e-1, 0.2, 1, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1, 10, 100}, ...
         'Ylim', [min(ylims(:,1))*0.8, max(ylims(:,2))*1.2], ....
         'fontsize', 15, ...
         'ticklength', [0.03, 0.03]);

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

function figure10(issave,folder_name,file_name)
% Figure 10: Spectra using different noise

params = utils.loadParameters_new;

noise_types = {'white', 'pink', 'brown', 'blue'};

init_x = 0.08; init_y = 0.15; spacing_x = 0.1; 
width = 0.4; height = 0.75;

h.fig = figure('Position', [200, 200, 800, 350]);

h.fmin = 0.01;
h.fmax = 1;
h.linewidth = 2;
h.linetypes = {'-', '--', '-.', ':'};
h.linecolors = [0 0 0; 1 0 1; 1 0 0; 0 0 1];
    
for noise_ind = 1:length(noise_types)
    out = utils.calc_spectra_noise(params, noise_types{noise_ind});
    f = out.w/(2*pi);
    fstart_ind = dsearchn(f', h.fmin);
    fend_ind = dsearchn(f', h.fmax);
    h.f = f(fstart_ind:fend_ind);
    normSpectra(noise_ind,:) = out.spectra(fstart_ind:fend_ind)/trapz(h.f, out.spectra(fstart_ind:fend_ind));
    noise_normSpectra(noise_ind,:) = out.noise_spectra(fstart_ind:fend_ind)/trapz(h.f, out.noise_spectra(fstart_ind:fend_ind)); 
end

% noise
h.ax = axes('Parent', h.fig, ...
             'Position', [init_x+(1-1)*spacing_x init_y width height]);
loglog(h.ax, NaN, NaN, 'HandleVisibility','off')
hold on
for noise_ind = 1:length(noise_types)
    loglog(h.ax, h.f, noise_normSpectra(noise_ind,:), h.linetypes{noise_ind}, 'linewidth', h.linewidth, ...
        'color', h.linecolors(noise_ind,:));
end
hold off 
set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 0.003, 1e-2, 0.05, 1e-1, 0.2, 1, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.003, 0.01, 0.05, 0.1, 0.2, 1, 10, 100}, ...
         'Ylim', [min(noise_normSpectra(:))*0.9, max(noise_normSpectra(:))*1.1], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]); 
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, '$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')
% title('noise spectrum', 'fontsize', 15) 

text(0.012, 0.6, 'white: $f^{0}$', 'interpreter', 'latex', 'fontsize', 18, ...
    'fontweight', 'b', 'color', h.linecolors(1,:))
text(0.11, 2.8, 'pink: $f^{-1}$', 'interpreter', 'latex', 'fontsize', 18, ...
    'fontweight', 'b', 'color', h.linecolors(2,:))
text(0.019, 29, 'brown: $f^{-2}$', 'interpreter', 'latex', 'fontsize', 18, ...
    'fontweight', 'b', 'color', h.linecolors(3,:))
text(0.031, 0.04, 'blue: $f^{1}$', 'interpreter', 'latex', 'fontsize', 18, ...
    'fontweight', 'b', 'color', h.linecolors(4,:))


% BOLD
h.ax = axes('Parent', h.fig, ...
             'Position', [init_x+width+(2-1)*spacing_x init_y width height]);
loglog(h.ax, NaN, NaN, 'HandleVisibility','off')
hold on
for noise_ind = 1:length(noise_types)
    loglog(h.ax, h.f, normSpectra(noise_ind,:), h.linetypes{noise_ind}, 'linewidth', h.linewidth, ...
        'color', h.linecolors(noise_ind,:));
end
hold off 
leg = legend(h.ax, 'white', 'pink', 'brown', 'blue');
set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 0.003, 1e-2, 0.05, 1e-1, 0.2, 1, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.003, 0.01, 0.05, 0.1, 0.2, 1, 10, 100}, ...
         'Ylim', [min(normSpectra(:))*0.9, max(normSpectra(:))*1.1], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);     
set(leg, 'fontSize', 15, 'location', 'southwest', 'box', 'off');
xlabel(h.ax, '$f$ (Hz)', 'fontsize', 15, 'interpreter', 'latex')
% ylabel(h.ax, '$P(f)$ (arbitrary units)', 'fontsize', 15, 'interpreter', 'latex')
% title('BOLD spectrum', 'fontsize', 15) 

annotation('textbox', [0.01 0.9, 0.1, 0.1], 'string', '(a)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'string', '(b)', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')

if issave
    set(h.fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(h.fig, '-painters', '-depsc', sprintf('%s/%s.eps',folder_name,file_name))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BASE FUNCTIONS TO PREPARE FIGURES AND DO CALCULATIONS
function prepare_figspectrum(h, normSpectra)

loglog(h.ax, h.f, normSpectra, h.linetype, 'linewidth', h.linewidth);
set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 0.003, 1e-2, 0.05, 1e-1, 0.2, 1, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.003, 0.01, 0.05, 0.1, 0.2, 1, 10, 100}, ...
         'Ylim', [min(normSpectra)*0.8, max(normSpectra)*1.2], ...
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);

function prepare_figVaryingSpectra(h, varying_param_normSpectra, what_spectra)

% spectra using nominal parameters
loglog(h.ax, NaN, NaN, 'HandleVisibility','off')
hold on;
% spectra using varying parameters within their ranges
for j = 1:h.number_points
    ylims(j, 1) = min(cellfun(@(x)min(varying_param_normSpectra{j}.(x)(:)), ...
                    fieldnames(varying_param_normSpectra{j})));
    ylims(j, 2) = max(cellfun(@(x)max(varying_param_normSpectra{j}.(x)(:)), ...
                    fieldnames(varying_param_normSpectra{j})));
    loglog(h.ax, h.f, varying_param_normSpectra{j}.(what_spectra), h.linetypes{j}, 'linewidth', h.linewidths(j), ...
           'color', h.linecolors(j,:));
end
hold off;

set(h.ax, 'Xlim', [h.fmin, h.fmax], ...
         'Xtick', [1e-5, 1e-4, 1e-3, 0.003, 1e-2, 0.05, 1e-1, 1, 5, 10, 100], ...
         'Xticklabel', {0.00001, 0.0001, 0.001, 0.003, 0.01, 0.05, 0.1, 1, 5, 10, 100}, ...
         'Ylim', [min(ylims(:,1))*0.8, max(ylims(:,2))*1.2], ....
         'fontsize', 15, ...
         'ticklength', [0.02, 0.02]);

function prepare_figBOLDspatiotemporal(h)

imagesc(h.ax, h.x*1e3, h.t, h.BOLD');
xlabel(h.ax, '$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
ylabel(h.ax, '$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')

colormap(colormap_bluetored)
colorbar('peer', h.ax)

set(h.ax, 'fontsize', 15, 'Xtick', [-20, -15, -10, -5, -2.5, 0, 2.5, 5, 10, 15, 20], ...
    'Xticklabel', {-20, -15, -10, -5, -2.5, 0 , 2.5, 5, 10, 15, 20}, ...
    'Ytick', [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], ...
    'Yticklabel', {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50}, ...
    'XLim', [min(h.x), max(h.x)*1.01]*1e3, 'YLim', [0, max(h.t)*1.01], ...
    'CLim', h.clims, 'Ydir', 'normal', 'ticklength', [0.02, 0.02])

function prepare_figBOLDtemporal(h)

plot(h.ax, h.t, h.BOLD_center, 'k-', 'linewidth', 1.5);

set(h.ax, 'Xtick', [0, 10, 20, 30, 40, 50], ...
    'Xticklabel', {0, 10, 20, 30, 40, 50},  ...
    'XLim', [0, max(h.t)*1.01], 'YLim', [min(h.BOLD_center)-0.02, max(h.BOLD_center)*1.2], ...
    'fontsize', 15, 'ticklength', [0.02, 0.02]); 
xlabel(h.ax, '$t$ (s)', 'fontsize', 15, 'interpreter', 'latex');

function h = calc_BOLD(h)
% calculate BOLD signal

d_max = 5e-3;
t_max = 20;

h.p.Nkx = 2^8;
h.p.Nw = 2^8;

h.x = d_max*(2*(0:h.p.Nkx-1)/h.p.Nkx - 1);
h.t = t_max*(2*(0:h.p.Nw-1)/h.p.Nw - 1);
h.center = dsearchn(h.x', 0);

kxsamp = (1/mean(diff(h.x)))*2*pi;
wsamp = (1/mean(diff(h.t)))*2*pi;

[kx, w] = generate_kw_1D(kxsamp, wsamp, h.p.Nkx, h.p.Nw);
T = calcTransFuncs_fromPhi_1D(kx, w, h.p);
neural_freq = 1;

h.BOLD = real(freq2coord_1D(T.T_Yphi.*neural_freq, kx, w));