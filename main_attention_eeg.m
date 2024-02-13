% Re-analysis of the XPAttention data at individual subject level. 
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

n_sub = 21; 

par.trial_dur = 33.6;   

par.roi_name = 'all'; % frontocentral, all
% par.roi_chans = {'F1', 'Fz', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'Cz', 'C2'}; 

par.ref_name = 'all'; % mastoid, all
% par.ref_chans = {'mast1', 'mast2'}; 


%% allocate table

col_names = {
    'subject', 'rhythm', 'task', ...
    'z_meter_fft_raw', ...
    'z_meter_fft_subtr', ...
    'z_meter_fft_sound', ...
    'z_snr' ...
    };

tbl_fft = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


col_names = {
    'subject', 'rhythm', 'task', ...
    'z_meter_acf_raw', ...
    'z_meter_acf_subtr', ...
    'z_meter_acf_sound' ...
    };

tbl_acf = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);



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


        %% process EEG

        % get acf
        % -------

        % withuout aperiodic subtraction    
        [acf, lags, ~, mX, freq] = get_acf(data, fs);    

        mX_subtracted = subtract_noise_bins(mX,...
                                        par.noise_bins(1), par.noise_bins(2)); 

        % with aperiodic subtraction    
        acf_subtracted = nan(size(acf)); 
        ap = nan(size(acf)); 
        
        parfor i_sub=1:size(data, 1)
            
            fprintf('getting 1/f-subtracted acf for sub-%02d\n', i_sub); 
            
            [acf_subtracted(i_sub, :, 1, 1, 1, :), ~, ...
             ap(i_sub, :, 1, 1, 1, :), ~, ~, ~, ~, ...
             optim_exitflag(i_sub, :)] = ...
                                get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'response_f0', par.response_f0, ...
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

        feat_acf = get_acf_features(acf, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 

        % get features for the raw spectra                                    
        feat_fft = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 

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
            repmat({task}, n_sub, 1), ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            repmat({feat_fft_coch.z_meter_rel}, n_sub, 1), ...
            num2cell(feat_fft.z_snr) ...
            ];
        

        tbl_fft = [tbl_fft; rows];
        
        rows = [...
            num2cell([1 : n_sub]'), ...
            repmat({rhythm}, n_sub, 1), ...
            repmat({task}, n_sub, 1), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            repmat({feat_acf_coch.z_meter_rel}, n_sub, 1)...
            ];

        tbl_acf = [tbl_acf; rows];
        
        
    end


end
    

%% save table

fname = sprintf('exp-attention_apFitMethod-%s_roi-%s_eegIndividual', ...
                par.ap_fit_method, par.roi_name); 
            
writetable(tbl_fft, fullfile(par.data_path, [fname, '_fft.csv'])); 
writetable(tbl_acf, fullfile(par.data_path, [fname, '_acf.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




