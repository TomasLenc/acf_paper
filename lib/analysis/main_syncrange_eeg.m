function main_syncrange_eeg(par)
% Re-analysis of the XPSyncRange data at individual subject level. 

rhythms = {'31', '26', '19', '37', '42', '6', '17', '41'}; 

n_rhythms = length(rhythms); 

%% allocate table

col_names = {
    'subject', 'rhythm', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% run

c = 1; 
data_to_plot = []; 

for i_rhythm=1:n_rhythms
    
    rhythm = rhythms{i_rhythm};
    
    fprintf('processing rhythm %s\n', rhythm);


    %% load data
    
    fpath_eeg = fullfile(par.eeg_path, 'syncrange'); 
    
    eeg = load(fullfile(fpath_eeg, ...
        sprintf('exp-syncrange_rhythm-%s_nTrials-9_eeg.mat', ...
                rhythm)));
    
    data = eeg.data;
    fs = eeg.fs;
    t = [0 : size(eeg.data, 2) - 1] / fs; 
        
    % make sure we don't have lags longer than half trial duration!
    trial_dur = size(data, 2) / fs; 

    par.lags_meter_rel = ...
        par.lags_meter_rel(par.lags_meter_rel < trial_dur/2); 

    par.lags_meter_unrel = ...
        par.lags_meter_unrel(par.lags_meter_unrel < trial_dur/2); 

    %% load stimulus
    
    fpath_stim = '/DATA1/XPSyncRange/ptb/eeg_7set/to_load'; 
    
    d = dir(fullfile(fpath_stim, sprintf('*rhythm%s_*.wav', rhythm)));
    
    [s, fs_s] = audioread(fullfile(fpath_stim, d.name));
    s = s(:, 1)'; 
    
    env = abs(hilbert(s)); 
    
    % low-pass filter for nicer acf plot 
    [b,a] = butter(2, 30/(fs_s/2), 'low'); 
    env = filtfilt(b, a, env); 

    %% process stimulus

    % get ACF (withuout aperiodic subtraction)
    [acf_coch, lags_coch, ~, mX_coch, freq_coch] = get_acf(env, fs_s);    

    % get ACF features
    feat_acf_coch = get_acf_features(...
                                acf_coch, lags_coch, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    

    % get features for the raw spectra                                    
    feat_fft_coch = get_fft_features(mX_coch, freq_coch, ...
                                     par.freq_meter_rel, par.freq_meter_unrel); 

    
    
    %% process EEG
    
    n_sub = size(data, 1); 
            
    % get acf
    % -------

    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(data, fs);    

    mX_subtracted = subtract_noise_bins(mX, ...
                            par.noise_bins(1),  par.noise_bins(2)); 

    % with aperiodic subtraction    
    acf_subtracted = nan(size(acf)); 

    parfor i_sub=1:size(data, 1)
        fprintf('sub-%02d\n', i_sub); 
        [acf_subtracted(i_sub, :)] = ...
                            get_acf(data(i_sub, :), fs, ...
                                   'rm_ap', true, ...
                                   'ap_fit_method', par.ap_fit_method, ...
                                   'response_f0', par.response_f0, ...
                                   'ap_fit_flims', par.ap_fit_flims, ...
                                   'verbose', false);
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


    % add features to table
    rows = [...
        num2cell([1 : n_sub]'), ...
        repmat({rhythm}, n_sub, 1), ...
        num2cell(feat_fft.z_meter_rel), ...
        num2cell(feat_acf.z_meter_rel), ...
        num2cell(feat_fft_subtracted.z_meter_rel), ...
        num2cell(feat_acf_subtracted.z_meter_rel), ...
        repmat({feat_fft_coch.z_meter_rel}, n_sub, 1), ...
        repmat({feat_acf_coch.z_meter_rel}, n_sub, 1), ...
        num2cell(feat_fft.z_snr) ...
        ];

    tbl = [tbl; rows];


    data_to_plot(c).rhythm = rhythm; 
    data_to_plot(c).mX = mX; 
    data_to_plot(c).mX_subtr = mX_subtracted;  
    data_to_plot(c).freq = freq; 
    data_to_plot(c).acf = acf; 
    data_to_plot(c).acf_subtr = acf_subtracted;
    data_to_plot(c).lags = lags; 
    data_to_plot(c).freq_coch = freq_coch; 
    data_to_plot(c).mX_coch = mX_coch; 
    data_to_plot(c).lags_coch = lags_coch; 
    data_to_plot(c).acf_coch = acf_coch; 
    c = c+1; 


end
    

fname = sprintf('exp-syncsweep_apFitMethod-%s_onlyHarm-%s_eegIndividual', ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics)); 

% save table
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save data 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 





