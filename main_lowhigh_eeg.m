function main_lowhigh_eeg(par)
% Re-analysis of the XPLowHigh data at individual subject level. 

load_path = fullfile(par.eeg_path, 'lowhigh'); 

rhythms = {'unsyncopated', 'syncopated'}; 
tones = {'L', 'H'}; 

n_rhythms = 2; 

rm_id = [1,4,6,11,14];

par.trial_dur = 50.4;  

par.roi_name = 'all'; 
% par.roi_chans = {'F1', 'Fz', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'Cz', 'C2'}; 

par.ref_name = 'all'; 


%% allocate table

col_names = {
    'subject', 'rhythm', 'tone', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);

data_to_plot = []; 

%% run

c = 1; 

for i_rhythm=1:n_rhythms
    
    for i_tone=1:2

        rhythm = rhythms{i_rhythm};
        tone = tones{i_tone}; 

        fprintf('processing rhythm: %s-%s\n', tone, rhythm);

        %% load data

        fname = fullfile(load_path, 'preprocessed', ...
            sprintf('avgRef_timeAvg_mergedParticipants %s_%s(19).lw6', ...
                    tone, rhythm));

        [header, data] = CLW_load(fname); 
        
        if strcmp(par.roi_name, 'all')
            par.roi_chans = {header.chanlocs.labels}; 
        end
        if strcmp(par.ref_name, 'all')
            par.ref_chans = {header.chanlocs.labels}; 
        end

        % remove bad subjects 
        [header, data] = RLW_arrange_epochs(header, data, ...
                                    setdiff([1:header.datasize(1)], rm_id)); 

        n_sub = header.datasize(1); 

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
                                            
        offset = 0; 

%         % if is unsyncopated, shift the ERP by 3 events (it's okay)
%         if strcmpi(rhythm_id, 'unsyncopated')
%             offset = 3 * 0.2; 
%         else
%             offset = 0; 
%         end

        [header, data] = RLW_segmentation(header, data, {'trig'}, ...
                            'x_start', 0+offset, 'x_duration', par.trial_dur); 

        % downsample 
        [header, data] = RLW_downsample(header, data, 'x_downsample_ratio', 16); 

        fs = 1/header.xstep; 
        t = [0 : header.datasize(end)-1] * header.xstep + header.xstart; 

        % make sure we don't have lags longer than half trial duration!
        par.lags_meter_rel = ...
            par.lags_meter_rel(par.lags_meter_rel < par.trial_dur/2); 

        par.lags_meter_unrel = ...
            par.lags_meter_unrel(par.lags_meter_unrel < par.trial_dur/2); 

        %% load stimulus

        fname = fullfile(load_path, 'Slaney_128coch_meddis_timeDomain_meanAcrossCF'); 

        coch_output = load(fname); % variables: freq, res_all, rowNames

        idx = ~cellfun(@isempty, ...
            strfind(coch_output.cond_names, [tone, '_standard_', rhythm])); 

        coch = coch_output.slaney(idx, :); 

        fs_coch = coch_output.fs;

        t_coch = coch_output.t; 

        % low-pass filter for nicer acf plot 
        [b,a] = butter(2, 30/(fs_coch/2), 'low'); 
        coch = filtfilt(b, a, coch); 

    %     % cut off the edges to remove filter artifacts 
    %     idx_start = round(2.4 * fs_coch); 
    %     N = round((par.trial_dur - 2 * 2.4) * fs_coch); 
    %     coch = coch(idx_start+1 : idx_start+N); 
    %     
    %     % ensure integer number of cycles 
    %     assert(mod(length(coch) / fs_coch, 2.4) == 0); 

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

        n_sub = size(data, 1); 

        % get acf
        % -------

        % withuout aperiodic subtraction    
        [acf, lags, ~, mX, freq] = get_acf(data, fs);    

        mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

        % with aperiodic subtraction    
        acf_subtracted = nan(size(acf)); 
        
        parfor i_sub=1:size(data, 1)
            
            fprintf('sub-%02d\n', i_sub); 
            
            [acf_subtracted(i_sub, :, 1, 1, 1, :)] = ...
                                get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'f0_to_ignore', par.f0_to_ignore, ...
                                       'ap_fit_flims', par.ap_fit_flims, ...
                                       'verbose', false);
        end

        % average across ROI channels 
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
            repmat({tone}, n_sub, 1), ...
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
        data_to_plot(c).tone = tone; 
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


end
    

%% save

fname = sprintf('exp-lowhigh_apFitMethod-%s_onlyHarm-%s_roi-%s_eegIndividual', ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics),...
                par.roi_name); 

% save table
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save data 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




