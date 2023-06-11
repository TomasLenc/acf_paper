function main_syncrange_eeg_boot(par)
% clear 
% par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%% parameters

n_boot = 1000;

n_trials_all = [9:-1:1]; 

rhythms = {'31', '26', '19', '37', '42', '6', '17', '41'}; 

n_rhythms = length(rhythms); 

% colors
cmap_name = 'Set1'; 
colors = num2cell(brewermap(n_rhythms, cmap_name), 2); 

%% allocate table

col_names = {
    'n_trials', 'rhythm', 'boot_sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl_boot = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


col_names = {
    'n_trials', 'rhythm', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl_grand = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% RUN

for i_n_trials=1:length(n_trials_all)
    
    n_trials = n_trials_all(i_n_trials); 

    for i_rhythm=1:n_rhythms
        
        rhythm_id = rhythms{i_rhythm};
        
        fprintf('processing rhythm %s\n', rhythm_id);
    
        %% load data
        
        fpath_eeg = fullfile(par.eeg_path, 'syncrange'); 
        
        eeg = load(fullfile(fpath_eeg, ...
            sprintf('exp-syncrange_rhythm-%s_nTrials-%d_eeg.mat', ...
                    rhythm_id, n_trials)));
        
        fs = eeg.fs;
        t = [0 : size(eeg.data, 2) - 1] / fs; 
                
        % make sure we don't have lags longer than half trial duration!
        trial_dur = size(eeg.data, 2) / fs; 
        
        par.lags_meter_rel = ...
            par.lags_meter_rel(par.lags_meter_rel < trial_dur/2); 
        
        par.lags_meter_unrel = ...
            par.lags_meter_unrel(par.lags_meter_unrel < trial_dur/2); 
        
        %% load stimulus
        
        fpath_stim = fullfile(par.eeg_path, 'syncrange', 'stimuli'); 
        
        d = dir(fullfile(fpath_stim, sprintf('*rhythm%s_*.wav', rhythm_id)));
        
        [s, fs_s] = audioread(fullfile(fpath_stim, d.name));
        s = s(:, 1)'; 
        env = abs(hilbert(s)); 
                
        t_s = [0 : length(s) - 1] / fs_s; 
                
        %% process stimulus
        
        % get ACF (withuout aperiodic subtraction)
        [acf_s, lags_s, ~, mX_s, freq_s] = get_acf(...
                                   env, fs_s);    
                               
        % get ACF features
        feat_acf_s = get_acf_features(acf_s, lags_s, ...
                               par.lags_meter_rel, par.lags_meter_unrel); 
                                   
        % get features for the raw spectra                      
        feat_fft_s = get_fft_features(mX_s, freq_s,...
                                   par.freq_meter_rel, par.freq_meter_unrel); 
        
        
        %% process bootstrap EEG
        
        n_sub = size(eeg.data, 1); 
        
        
        % let's do the boostarp in batches
        batch_size = 500; 
        
        n_boot_batches = ceil(n_boot / batch_size); 
        
        % allocate 
        feat_acf_boot = struct('z_meter_rel', []); 
        feat_acf_subtracted_boot = struct('z_meter_rel', []); 
        feat_fft_boot = struct('z_meter_rel', [], 'z_snr', []); 
        feat_fft_subtracted_boot = struct('z_meter_rel', []); 
        feat_ap_boot = struct('offset', [], 'exponent', []); 
        
        for i_batch=1:n_boot_batches
            
            n_boot_batch = min(n_boot - (i_batch-1)*batch_size, batch_size);
            
            eeg_boot = nan(n_boot_batch, size(eeg.data, 2));

            for i_boot=1:n_boot_batch
                idx = randsample(n_sub, n_sub, true);
                eeg_boot(i_boot, :) = mean(eeg.data(idx, :), 1); 
            end

            % get acf
            % -------

            % withuout aperiodic subtraction    
            [acf, lags, ~, mX, freq] = get_acf(eeg_boot, fs);    

            mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

            % with aperiodic subtraction    
            [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_exitflag] = ...
                                        get_acf(eeg_boot, fs, ...
                                               'rm_ap', true, ...
                                               'f0_to_ignore', 1 / 2.4, ...
                                               'min_freq', 0.1, ...
                                               'max_freq', 9);      
            if any(~optim_exitflag)
                warning('ap-fit didnt converge %d/%d reps', sum(~optim_exitflag), n_rhythms); 
            end
            
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
            
            % 1/f parameter estimates
            feat_ap_boot.offset = cat(1, ...
                feat_ap_boot.offset, cellfun(@(x) x(1), par_ap));  
            
            feat_ap_boot.exponent = cat(1, ...
                feat_ap_boot.exponent, cellfun(@(x) x(2), par_ap));            


        end
        
        
        %% process EEG grand average 
        
        % take grand average across subjects
        eeg_grand = mean(eeg.data, 1); 
        
        % get acf
        % -------
                                           
        % withuout aperiodic subtraction    
        [acf_grand, lags, ~, mX_grand, freq] = get_acf(eeg_grand, fs);    
                                               
        % with aperiodic subtraction    
        [acf_subtracted_grand, ~, ap_grand, ~, ~, par_ap] = ...
                                    get_acf(eeg_grand, fs, ...
                                           'rm_ap', true, ...
                                           'f0_to_ignore', 1 / 2.4, ...
                                           'min_freq', 0.1, ...
                                           'max_freq', 9);
        
        feat_ap.offset = cellfun(@(x) x(1), par_ap);            
        feat_ap.exponent = cellfun(@(x) x(2), par_ap);            
                                       
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
    
        mX_subtracted_grand = subtract_noise_bins(...
                                mX_grand, par.noise_bins(1), par.noise_bins(2)); 

        feat_fft_subtracted = get_fft_features(mX_subtracted_grand, freq, ...
                               par.freq_meter_rel, par.freq_meter_unrel); 
                                   
                                           
        %% plot example 
        
        if i_n_trials == 1
            if i_rhythm==1
                f = figure('color','white', ...
                           'position', [95, 67, 1062, 170 * n_rhythms]); 
                pnl_example = panel(f); 
                pnl_example.pack('v', n_rhythms); 
                pnl_example.margin = [5, 10, 25, 25]; 
            end
            plot_example(eeg_grand, t, ...
                             acf_grand, lags, ...
                             ap_grand, ...
                             mX_grand, freq, ...
                             par.lags_meter_rel, par.lags_meter_unrel, ...
                             par.freq_meter_rel, par.freq_meter_unrel, ...
                             'pnl', pnl_example(i_rhythm), ...
                             'subplot_proportions', [35, 10, 55], ...
                             'min_lag', 0.2, ...
                             'max_lag', par.max_lag, ...
                             'max_freq', par.max_freq_plot, ...
                             'plot_time_xaxis', i_rhythm == n_rhythms, ...
                             'plot_xlabels', i_rhythm == n_rhythms, ...
                             'plot_xticks', i_rhythm == n_rhythms, ...
                             'plot_features', false, ...
                             'mX_subtr', mX_subtracted_grand, ...
                             'acf_subtr', acf_subtracted_grand, ...
                             'time_col', colors{i_rhythm}, ...
                             'prec', 1e6, ...
                             'fontsize', par.fontsize, ...
                             'normalize_acf_for_plotting', false);                                        
            f.Name = rhythms{i_rhythm};     
            pnl_example(i_rhythm).margintop = 35; 
            
        end
    
        %% add features to table  
        
        rows = [...
            repmat({n_trials}, n_boot, 1), ...
            repmat({rhythm_id}, n_boot, 1), ...
            num2cell([1:n_boot]'), ...
            num2cell(feat_fft_boot.z_meter_rel), ...
            num2cell(feat_acf_boot.z_meter_rel), ...
            num2cell(feat_fft_subtracted_boot.z_meter_rel), ...
            num2cell(feat_acf_subtracted_boot.z_meter_rel), ...
            repmat({feat_fft_s.z_meter_rel}, n_boot, 1), ...
            repmat({feat_acf_s.z_meter_rel}, n_boot, 1), ...
            num2cell(feat_fft_boot.z_snr), ...
            num2cell(feat_ap_boot.offset), ...
            num2cell(feat_ap_boot.exponent) ...
            ];
        
        tbl_boot = [tbl_boot; rows];
        
        rows = [...
            {n_trials}, ...
            {rhythm_id}, ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            {feat_fft_s.z_meter_rel}, ...
            {feat_acf_s.z_meter_rel}, ...
            num2cell(feat_fft.z_snr), ...
            num2cell(feat_ap.offset), ...
            num2cell(feat_ap.exponent) ...
            ];
        
        tbl_grand = [tbl_grand; rows];      
        
    end
        
    if i_n_trials == 1
        fname = sprintf('exp-syncrange_response-eeg_nTrials-%d_boot_examples.svg', n_trials); 
        save_fig(f, fullfile(par.data_path, fname))
    end
    
end

%% save tables

fname = sprintf('exp-syncrange_eegGrand'); 
writetable(tbl_grand, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 

fname = sprintf('exp-syncrange_eegBoot'); 
writetable(tbl_boot, fullfile(par.data_path, [fname, '.csv'])); 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 









