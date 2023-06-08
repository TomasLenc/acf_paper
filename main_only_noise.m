function main_only_noise()
% clear

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%% simulate

noise_exponent = -1.5; 

fit_knee = false; 

noise_type = 'eeg'; % eeg, fractal

% number of simulated repetitions 
n_rep = 500; % 1000

%% generate noise

N = round(par.trial_dur * par.fs); 

% generate noisy signal (simulataneously for all repetitions) 
if strcmp(noise_type, 'fractal')

    x = get_colored_noise2([n_rep, N], par.fs, noise_exponent); 

elseif strcmp(noise_type, 'eeg')

    x = prepare_eeg_noise(n_rep, par.trial_dur); 

else
    error('noise type "%s" not implemented', noise_type);
end


%% analyse

% get acf
% -------

% withuout aperiodic subtraction    
[acf, lags, ~, mX, freq] = get_acf(x, par.fs);    

mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

% with aperiodic subtraction    
[acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr] = get_acf(x, par.fs, ...
                                   'rm_ap', true, ...
                                   'f0_to_ignore', 1/2.4, ...
                                   'min_freq', 0.1, ...
                                   'max_freq', 9);  

feat_ap.offset = cellfun(@(x) x(1), par_ap);            
feat_ap.exponent = cellfun(@(x) x(2), par_ap);            

% get features
% ------------

feat_acf = get_acf_features(acf, lags, ...
                            par.lags_meter_rel, par.lags_meter_unrel);    

feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                             par.lags_meter_rel, par.lags_meter_unrel); 

% get features for the raw spectra                                    
tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
feat_fft.z_meter_rel = tmp.z_meter_rel; 

feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                   par.noise_bins_snr(1), ...
                                   par.noise_bins_snr(2)); 

% get features for the 1/f-subtracted spectra                                    
feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                       par.freq_meter_rel, par.freq_meter_unrel);




%% save 

if strcmp(noise_type, 'fractal')
    fname = sprintf('onlyNoise_noise-fractal_exp-%.1f_nrep-%d', ...
                   noise_exponent, n_rep); 
elseif strcmp(noise_type, 'eeg')
    fname = sprintf('onlyNoise_noise-eeg_nrep-%d', ...
                   n_rep); 
else
    error('noise type "%s" not implemented', noise_type);
end

tbl = [
    num2cell(feat_fft.z_meter_rel), ...
    num2cell(feat_fft_subtracted.z_meter_rel), ...
    num2cell(feat_acf.z_meter_rel), ...
    num2cell(feat_acf_subtracted.z_meter_rel)
    ];

tbl = cell2table(tbl, ...
    'VariableNames', {'fft_raw', 'fft_subtr', 'acf_raw', 'acf_subtr'}); 

writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 

























