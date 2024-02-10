function main_lowhigh_tap(par)
% Re-analysis of the XPLowHigh data at individual subject level. 

load_path = fullfile(par.eeg_path, 'lowhigh'); 

rhythms = {'unsyncopated', 'syncopated'}; 
tones = {'L', 'H'}; 

n_rhythms = 2; 

rm_id = [1,4,6,11,14];

keep_id = setdiff([1:19], rm_id); 

par.trial_dur = 57.6;  

N = round(44100 * par.trial_dur); 


%% allocate table

col_names = {
    'subject', 'trial', 'rhythm', 'tone', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);

data_to_plot = []; 

%% run

c = 1; 

for i_sub=1:length(keep_id)
    
    sub = keep_id(i_sub); 

    % load data
    d = dir(fullfile(load_path, 'tapping', ...
        sprintf('%d_*.mat', sub))); 

    load(fullfile(d.folder, d.name)); 

    assert(strcmp(res.subjectID, num2str(sub)));

    for i_rhythm=1:n_rhythms

        for i_tone=1:2

            rhythm = rhythms{i_rhythm};
            tone = tones{i_tone}; 

            fprintf('processing sub-%02d rhythm: %s-%s\n', sub, tone, rhythm);

            cond_mask = ~cellfun(@isempty, ...
                strfind(res.blockOrder, sprintf('%s_%s', tone, rhythm))); 
            
            data1 = res.data.block(cond_mask).trial(1).tapData(1, :); 
            data2 = res.data.block(cond_mask).trial(2).tapData(1, :); 
            
            % filter
            [b,a] = butter(2, 30/(44100/2), 'low'); 
            data1 = filtfilt(b, a, data1); 
            data2 = filtfilt(b, a, data2); 
            
            % merge
            data = [data1(1:N); data2(1:N)]; 

            % downsample by factor of 100
            data = data(:, 1 : 100 : end); 

            fs = 44100 / 100; 
            t = [0 : size(data, 2) - 1] / fs; 

            % make sure we don't have lags longer than half trial duration!
            par.lags_meter_rel = ...
                par.lags_meter_rel(par.lags_meter_rel < par.trial_dur/2); 

            par.lags_meter_unrel = ...
                par.lags_meter_unrel(par.lags_meter_unrel < par.trial_dur/2); 

            %% process tapping

            n_trials = size(data, 1); 

            % get acf
            % -------

            % withuout aperiodic subtraction    
            [acf, lags, ~, mX, freq] = get_acf(data, fs);    

            mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

            % with aperiodic subtraction    
            acf_subtracted = nan(size(acf)); 

            parfor i_trial=1:size(data, 1)
                [acf_subtracted(i_trial, :)] = ...
                                    get_acf(data(i_trial, :), fs, ...
                                           'rm_ap', true, ...
                                           'ap_fit_method', par.ap_fit_method, ...
                                           'f0_to_ignore', par.f0_to_ignore, ...
                                           'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                                           'keep_band_around_f0_harmonics', par.ap_band_around_harmonics, ...
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
                repmat({sub}, n_trials, 1), ...
                num2cell([1 : n_trials]'), ...
                repmat({rhythm}, n_trials, 1), ...
                repmat({tone}, n_trials, 1), ...
                num2cell(feat_fft.z_meter_rel), ...
                num2cell(feat_acf.z_meter_rel), ...
                num2cell(feat_fft_subtracted.z_meter_rel), ...
                num2cell(feat_acf_subtracted.z_meter_rel), ...
                num2cell(feat_fft.z_snr) ...
                ];

            tbl = [tbl; rows];


            data_to_plot(c).sub = sub; 
            data_to_plot(c).rhythm = rhythm; 
            data_to_plot(c).tone = tone; 
            data_to_plot(c).mX = mX; 
            data_to_plot(c).mX_subtr = mX_subtracted;  
            data_to_plot(c).freq = freq; 
            data_to_plot(c).acf = acf; 
            data_to_plot(c).acf_subtr = acf_subtracted;
            data_to_plot(c).lags = lags; 

            c = c+1; 

        end

    end
    
end


%% save

fname = sprintf('exp-lowhigh_apFitMethod-%s_onlyHarm-%s_tapIndividual', ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics)); 

% save table
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save data 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




