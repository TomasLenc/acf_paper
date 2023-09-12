function main_only_noise(par, varargin)
% clear
% par = get_par(); 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); 

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;

%% generate noise


if size(noise, 1) < par.n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          par.n_rep, size(noise, 1)); 
else
    noise = noise(1:par.n_rep, :); 
end

%% analyse

% get acf
% -------

% withuout aperiodic subtraction    
[acf, lags, ~, mX, freq] = get_acf(noise, par.fs);    

mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

% with aperiodic subtraction    
acf_subtracted = nan(size(acf)); 

parfor i_rep=1:par.n_rep
    fprintf('analysing %d/%d\n', i_rep, par.n_rep); 
    
    [acf_subtracted(i_rep, :)] = get_acf(noise(i_rep, :), par.fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'f0_to_ignore', par.f0_to_ignore, ...
                                       'ap_fit_flims', par.ap_fit_flims ...
                                       ); 
end

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

fname = sprintf('ir-%s_noise-%s_apFitMethod-%s_onlyHarm-%s_onlyNoise', ...
               par.ir_type, ...
               par.noise_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 

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

























