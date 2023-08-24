% function main_syncrange_eeg(par)
clear 
par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))


%% parameters

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.05; 

load_path = '/datadisk/projects/XPLowHigh'; 

rhythms = {'unsyncopated', 'syncopated'}; 
tones = {'L', 'H'}; 

n_rhythms = 2; 

rm_id = [1, 4, 6, 11, 14];

trial_dur = 50.4;  

n_boot = 200;


%% allocate table

col_names = {
    'rhythm', 'tone', 'boot_sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl_boot = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);

col_names = {
    'rhythm', 'tone', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl_grand = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% lags of interest 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
min_lag = 0;
max_lag = trial_dur / 2; 

lag_base_incl_meter_rel = [0.8]; 
lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
lag_base_excl_meter_unrel = [0.8]; 

par.lags_meter_rel = get_lag_harmonics(...
                            lag_base_incl_meter_rel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_rel ...
                            ); 
                        
par.lags_meter_unrel = get_lag_harmonics(...
                            lag_base_incl_meter_unrel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_unrel ...
                            ); 

%% run

for i_rhythm=1:n_rhythms
    
    for i_tone=1:2
        
        rhythm_id = rhythms{i_rhythm};
        tone_id = tones{i_tone}; 

        fprintf('processing rhythm: %s-%s\n', tone_id, rhythm_id);


        %% load data

        fname = fullfile(load_path, 'preprocessed', ...
            sprintf('avgRef_timeAvg_mergedParticipants %s_%s(19).lw6', ...
                    tone_id, rhythm_id));

        [header, data] = CLW_load(fname); 

        % remove bad subjects 
        [header, data] = RLW_arrange_epochs(header, data, ...
                                    setdiff([1:header.datasize(1)], rm_id)); 

        n_sub = header.datasize(1); 

        n_chan = header.datasize(2); 

        % filter
        [header, data] = RLW_butterworth_filter(header, data, ...
                                                'filter_type', 'lowpass', ...
                                                'high_cutoff', 30, ...
                                                'filter_order', 2); 

        % if is unsyncopated, shift the ERP by 3 events!!!
        if strcmpi(rhythm_id, 'unsyncopated')
            offset = 3 * 0.2; 
        else
            offset = 0; 
        end

        [header, data] = RLW_segmentation(header, data, {'trig'}, ...
                            'x_start', 0+offset, 'x_duration', trial_dur); 

        % downsample 
        [header, data] = RLW_downsample(header, data, 'x_downsample_ratio', 16); 

        fs = 1/header.xstep; 
        t = [0 : header.datasize(end)-1] * header.xstep + header.xstart; 

        % make sure we don't have lags longer than half trial duration!
        par.lags_meter_rel = ...
            par.lags_meter_rel(par.lags_meter_rel < trial_dur/2); 

        par.lags_meter_unrel = ...
            par.lags_meter_unrel(par.lags_meter_unrel < trial_dur/2); 

        %% load stimulus

        fname = fullfile(load_path, 'Slaney_128coch_meddis_timeDomain_meanAcrossCF'); 

        coch_output = load(fname); % variables: freq, res_all, rowNames

        idx = ~cellfun(@isempty, ...
            strfind(coch_output.cond_names, [tone_id, '_standard_', rhythm_id])); 

        coch = coch_output.slaney(idx, :); 

        fs_coch = coch_output.fs;

        t_coch = coch_output.t; 

        % low-pass filter for nicer acf plot 
        [b,a] = butter(2, 30/(fs_coch/2), 'low'); 
        coch = filtfilt(b, a, coch); 

        % ensure eeg and coch have the same duration 
        assert(length(coch) / fs_coch - length(data) / fs < 1e-4)


        %% process stimulus

        % get ACF (withuout aperiodic subtraction)
        [acf_coch, lags_coch, ~, mX_coch, freq_coch] = get_acf(coch, fs_coch);    

        % get ACF features
        feat_acf_coch = get_acf_features(...
                                    acf_coch, lags_coch, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        % get features for the raw spectra                                    
        feat_fft_coch = get_fft_features(mX_coch, freq_coch, ...
                                         par.freq_meter_rel, par.freq_meter_unrel); 


        %% process EEG 
        
        % let's do the boostarp in batches
        batch_size = 50; 
        
        n_boot_batches = ceil(n_boot / batch_size); 
        
        % allocate 
        feat_acf_boot = struct('z_meter_rel', []); 
        feat_acf_subtracted_boot = struct('z_meter_rel', []); 
        feat_fft_boot = struct('z_meter_rel', [], 'z_snr', []); 
        feat_fft_subtracted_boot = struct('z_meter_rel', []); 
        
        for i_batch=1:n_boot_batches
                        
            n_boot_batch = min(n_boot - (i_batch-1)*batch_size, batch_size);
            
            fprintf('processing batch %d/%d (%d samples)\n', ...
                    i_batch, n_boot_batches, n_boot_batch); 

            eeg_boot = nan(n_boot_batch, size(data, 2), size(data, 6));

            for i_boot=1:n_boot_batch
                
                idx = randsample(n_sub, n_sub, true);
                
                eeg_boot(i_boot, :, :) = mean(data(idx, :, 1, 1, 1, :), 1); 
                
            end

            % get acf
            % -------

            % withuout aperiodic subtraction    
            [acf, lags, ~, mX, freq] = get_acf(eeg_boot, fs);    

            mX_subtracted = subtract_noise_bins(mX, ...
                                    par.noise_bins(1),  par.noise_bins(2)); 

            % with aperiodic subtraction    
            [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_exitflag] = ...
                                        get_acf(eeg_boot, fs, ...
                                               'rm_ap', true, ...
                                               'f0_to_ignore', 1 / 2.4, ...
                                               'ap_fit_flims', [0.1, 9], ...
                                               'plot_diagnostic', false);      
            if any(~optim_exitflag)
                warning('ap-fit didnt converge %d/%d reps', sum(~optim_exitflag), n_rhythms); 
            end
            
            % average across channels 
            mX = squeeze(mean(mX, 2)); 
            mX_subtracted = squeeze(mean(mX_subtracted, 2)); 
            acf = squeeze(mean(acf, 2)); 
            acf_subtracted = squeeze(mean(acf_subtracted, 2)); 

            % get features
            % ------------

            % raw ACF                          
            tmp = get_acf_features(acf, lags, ...
                                   par.lags_meter_rel, par.lags_meter_unrel); 
                                        
            feat_acf_boot.z_meter_rel = cat(1, ...
                            feat_acf_boot.z_meter_rel, ...
                            tmp.z_meter_rel); 

            % noise-subtracted ACF                            
            tmp = get_acf_features(acf_subtracted, lags, ...
                                         par.lags_meter_rel, par.lags_meter_unrel); 
                                     
            feat_acf_subtracted_boot.z_meter_rel = cat(1, ...
                            feat_acf_subtracted_boot.z_meter_rel, ...
                            tmp.z_meter_rel); 

                        
            % raw FFT              
            tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
            
            feat_fft_boot.z_meter_rel = cat(1, ...
                feat_fft_boot.z_meter_rel, tmp.z_meter_rel); 

            tmp = get_z_snr(mX, freq, par.frex, ...
                            par.noise_bins_snr(1), par.noise_bins_snr(2)); 
            
            feat_fft_boot.z_snr = cat(1, feat_fft_boot.z_snr, tmp); 

            % noise-subtracted FFT
            tmp = get_fft_features(mX_subtracted, freq, ...
                                   par.freq_meter_rel, par.freq_meter_unrel);
            
            feat_fft_subtracted_boot.z_meter_rel = cat(1, ...
                feat_fft_subtracted_boot.z_meter_rel, tmp.z_meter_rel); 

        end
        
        
        %% process EEG grand average 
        
        % take grand average across subjects
        data_grand = mean(data, 1); 
        
        % get acf
        % -------
                                           
        % withuout aperiodic subtraction    
        [acf_grand, lags, ~, mX_grand, freq] = get_acf(data_grand, fs);    
                                               
        % with aperiodic subtraction    
        [acf_subtracted_grand, ~, ap_grand, ~, ~, par_ap] = ...
                                    get_acf(data_grand, fs, ...
                                           'rm_ap', true, ...
                                           'f0_to_ignore', 1 / 2.4, ...
                                           'ap_fit_flims', [0.1, 9], ...
                                           'plot_diagnostic', false);
                                        
        mX_subtracted_grand = subtract_noise_bins(...
                                mX_grand, par.noise_bins(1), par.noise_bins(2)); 

        % average across channels 
        mX_grand = ensure_row(squeeze(mean(mX_grand, 2))); 
        mX_subtracted_grand = ensure_row(squeeze(mean(mX_subtracted_grand, 2))); 
        acf_grand = ensure_row(squeeze(mean(acf_grand, 2))); 
        acf_subtracted_grand = ensure_row(squeeze(mean(acf_subtracted_grand, 2))); 

        % get features
        % ------------

        % ACF features                      
        feat_acf = get_acf_features(acf_grand, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel); 
                                    
        feat_acf_subtracted = get_acf_features(acf_subtracted_grand, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 
                                 
                                     
        % FFT features         
        feat_fft = get_fft_features(mX_grand, freq, ...
                              par.freq_meter_rel, par.freq_meter_unrel); 
                          
        feat_fft.z_snr = get_z_snr(mX_grand, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 
    
        feat_fft_subtracted = get_fft_features(mX_subtracted_grand, freq, ...
                               par.freq_meter_rel, par.freq_meter_unrel); 
        
                           
        %% add features to table  
        
        rows = [...
            repmat({rhythm_id}, n_boot, 1), ...
            repmat({tone_id}, n_boot, 1), ...
            num2cell([1:n_boot]'), ...
            num2cell(feat_fft_boot.z_meter_rel), ...
            num2cell(feat_acf_boot.z_meter_rel), ...
            num2cell(feat_fft_subtracted_boot.z_meter_rel), ...
            num2cell(feat_acf_subtracted_boot.z_meter_rel), ...
            repmat({feat_fft_coch.z_meter_rel}, n_boot, 1), ...
            repmat({feat_acf_coch.z_meter_rel}, n_boot, 1), ...
            num2cell(feat_fft_boot.z_snr) ...
            ];
        
        tbl_boot = [tbl_boot; rows];
        
        rows = [...
            {rhythm_id}, ...
            {tone_id}, ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            {feat_fft_coch.z_meter_rel}, ...
            {feat_acf_coch.z_meter_rel}, ...
            num2cell(feat_fft.z_snr) ...
            ];
        
        tbl_grand = [tbl_grand; rows];      
                                                      
                           
    end


end
    
%% save tables

fname = sprintf('exp-lowhigh_eegGrand'); 
writetable(tbl_grand, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 

fname = sprintf('exp-lowhigh_eegBoot'); 
writetable(tbl_boot, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 





