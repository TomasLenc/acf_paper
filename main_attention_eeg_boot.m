% function main_syncrange_eeg(par)
clear 
par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))


%% parameters

load_path = '/datadisk/projects/Attention'; 

rhythms = {'unsyncopated', 'syncopated'}; 
tasks = {'tempo', 'pitch', 'arithmetics'}; 

n_rhythms = length(rhythms); 
n_tasks = length(tasks); 

par.trial_dur = 33.6;   

par.roi_name = 'all'; % frontocentral, all
% par.roi_chans = {'F1', 'Fz', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'Cz', 'C2'}; 

par.ref_name = 'all'; % mastoid, all
par.ref_chans = {'mast1', 'mast2'}; 

par.n_boot = 500;

%% allocate table

col_names = {
    'rhythm', 'task', 'boot_sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl_boot = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


col_names = {
    'rhythm', 'task', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl_grand = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);



%% frequencies of interest 


%% lags of interest 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
par.max_lag = par.trial_dur / 2; 

par.lag_base_incl_meter_rel = [0.8]; 
par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
par.lag_base_excl_meter_unrel = [0.4]; 

par.lags_meter_rel = get_lag_harmonics(...
                            par.lag_base_incl_meter_rel, ...
                            par.max_lag,...
                            'lag_harm_to_exclude', par.lag_base_excl_meter_rel ...
                            ); 
                        
par.lags_meter_unrel = get_lag_harmonics(...
                            par.lag_base_incl_meter_unrel, ...
                            par.max_lag,...
                            'lag_harm_to_exclude', par.lag_base_excl_meter_unrel ...
                            ); 

%% run

for i_rhythm=1:n_rhythms
    
    for i_task=1:n_tasks

        rhythm = rhythms{i_rhythm};
        task = tasks{i_task}; 

        fprintf('processing rhythm: %s-%s\n', rhythm, task);

        %% get EEG data
        
        fname = fullfile(load_path, 'preprocessed', ...
            sprintf('%s-%s timeAvg mergedParticipants*.lw6', ...
                    rhythm, task));
                
        d = dir(fullfile(fname)); 

        [header, data] = CLW_load(fullfile(d.folder, d.name)); 
        
        [header, data] = RLW_arrange_epochs(header, data, ...
                                            [5 : header.datasize(1)]); 
        
        n_sub = header.datasize(1); 

        n_chan = header.datasize(2); 

        if strcmp(par.roi_name, 'all')
            par.roi_chans = {header.chanlocs.labels}; 
        end
        if strcmp(par.ref_name, 'all')
            par.ref_chans = {header.chanlocs.labels}; 
        end
        
        % reference
        [header, data] = RLW_rereference(header, data,...
                                'apply_list', {header.chanlocs.labels}, ...
                                'reference_list', par.ref_chans); 
                            
        % select channels of interest 
        [header, data] = RLW_arrange_channels(header, data, par.roi_chans); 
        
        assert(header.datasize(2) == length(par.roi_chans)); 
        
        % average channels of interest (unless we're taking all channels)
        if ~strcmp(par.roi_name, 'all')
            [header, data] = RLW_pool_channels(header, data, par.roi_chans, ...
                                               'keep_original_channels', 0); 
        end
        
        % filter
        [header, data] = RLW_butterworth_filter(header, data, ...
                                                'filter_type', 'lowpass', ...
                                                'high_cutoff', 30, ...
                                                'filter_order', 2); 

        % segment                                     
        trig = unique({header.events.code}); 
        assert(length(trig) == 1); 
        
        [header, data] = RLW_segmentation(header, data, trig, ...
                            'x_start', 0, 'x_duration', par.trial_dur); 

        % downsample 
        [header, data] = RLW_downsample(header, data, 'x_downsample_ratio', 4); 

        fs = 1/header.xstep; 
        t = [0 : header.datasize(end)-1] * header.xstep + header.xstart; 

        % make sure we don't have lags longer than half trial duration!
        par.lags_meter_rel = ...
            par.lags_meter_rel(par.lags_meter_rel < par.trial_dur/2); 

        par.lags_meter_unrel = ...
            par.lags_meter_unrel(par.lags_meter_unrel < par.trial_dur/2); 


        %% load stimulus

        d = dir(fullfile(load_path, 'urear', sprintf('UREAR_*_%s.mat', rhythm))); 
        
        coch_data = load(fullfile(d.folder, d.name)); 
        
        coch = sum(coch_data.output.AN.an_sout, 1); 
        
        fs_coch = coch_data.output.AN.fs; 
        
        t_coch = [0 : length(coch)-1] / fs_coch; 

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


        %% bootstrap EEG 
        
        % let's do the boostarp in batches
        batch_size = 50; 
        
        n_boot_batches = ceil(par.n_boot / batch_size); 
        
        % allocate 
        feat_acf_boot = struct('z_meter_rel', []); 
        feat_acf_subtracted_boot = struct('z_meter_rel', []); 
        feat_fft_boot = struct('z_meter_rel', [], 'z_snr', []); 
        feat_fft_subtracted_boot = struct('z_meter_rel', []); 
        
        for i_batch=1:n_boot_batches
                        
            n_boot_batch = min(par.n_boot - (i_batch-1)*batch_size, batch_size);
            
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
            acf_subtracted = nan(size(acf)); 
            ap = nan(size(acf)); 
            
            parfor i_boot=1:size(eeg_boot, 1)
                
                fprintf('boot-%02d\n', i_boot); 
                
                [acf_subtracted(i_boot, :, :), ~, ...
                 ap(i_boot, :, :), ~, ~, ~, ~, ...
                 optim_exitflag(i_boot, :)] = ...
                                    get_acf(eeg_boot(i_boot, :, :), fs, ...
                                           'rm_ap', true, ...
                                           'ap_fit_method', par.ap_fit_method, ...
                                           'f0_to_ignore', par.f0_to_ignore, ...
                                           'ap_fit_flims', par.ap_fit_flims, ...
                                           'verbose', false);
            end
            
            if strcmp(par.ap_fit_method, 'fooof') && any(~optim_exitflag)
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
                                           'ap_fit_method', par.ap_fit_method, ...
                                           'f0_to_ignore', par.f0_to_ignore, ...
                                           'ap_fit_flims', par.ap_fit_flims, ...
                                           'plot_diagnostic', false, ...
                                           'verbose', true);
                                        
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
            repmat({rhythm}, par.n_boot, 1), ...
            repmat({task}, par.n_boot, 1), ...
            num2cell([1:par.n_boot]'), ...
            num2cell(feat_fft_boot.z_meter_rel), ...
            num2cell(feat_acf_boot.z_meter_rel), ...
            num2cell(feat_fft_subtracted_boot.z_meter_rel), ...
            num2cell(feat_acf_subtracted_boot.z_meter_rel), ...
            repmat({feat_fft_coch.z_meter_rel}, par.n_boot, 1), ...
            repmat({feat_acf_coch.z_meter_rel}, par.n_boot, 1), ...
            num2cell(feat_fft_boot.z_snr) ...
            ];
        
        tbl_boot = [tbl_boot; rows];
        
        rows = [...
            {rhythm}, ...
            {task}, ...
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

%% save table

fname = sprintf('exp-attention_apFitMethod-%s_roi-%s_eegGrand', ...
                par.ap_fit_method, par.roi_name); 
writetable(tbl_grand, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 

fname = sprintf('exp-attention_apFitMethod-%s_roi-%s_eegBoot', ...
                par.ap_fit_method, par.roi_name); 
writetable(tbl_boot, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




