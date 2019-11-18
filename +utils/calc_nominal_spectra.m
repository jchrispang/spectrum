function [nominal_spectra, nominal_normSpectra] = calc_nominal_spectra(f)

params = utils.loadParameters_new;
[nominal_spectra, nominal_normSpectra] = utils.calc_spectra(f, params);