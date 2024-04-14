function main_infant_eeg(par)
% Re-analysis of the XPInfant data at individual subject level. 

%% parameters

load_path = fullfile(par.eeg_path, 'infant'); 

rhythms = {'unsync', 'sync'}; 
tones = {'low', 'high'}; 

n_rhythms = length(rhythms); 

par.trial_dur = 60;   

subjects = [1:20]; 

% mastoids 
ref_chans = {'E57','E100'}; 

% set this to true if you want the preprocessed trials to be chunked before
% averaging (this might help to increase the SNR?)
preproc_chunk_do = true; 
preproc_chunk_onset = 2.4; 
preproc_chunk_dur = 26.4; 

n_chan = 92; 

front_chan = {...
        'E46', 'E41', 'E36',  'E40', 'E35',  'E29',  'E34', 'E28', 'E20', 'E27', 'E24', 'E19', ...
        'E102', 'E103', 'E104', 'E109','E110','E111', 'E116','E117', 'E118', 'E123', 'E124', 'E4'}; 

roi = 'front'; 

%% allocate table

col_names = {
    'subject', 'rhythm', 'tone', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% frequencies of interest 

% exclude first two harmonics as in the original paper
par.freq_meter_rel = par.freq_meter_rel(par.freq_meter_rel > 1/2.4); 
par.freq_meter_unrel = par.freq_meter_unrel(par.freq_meter_unrel > 1/2.4); 
par.frex = par.frex(par.frex > 1/2.4); 

%% run

c = 1; 

for i_rhythm=1:n_rhythms
    
    for i_tone=1:2

        rhythm = rhythms{i_rhythm};
        tone = tones{i_tone}; 

        fprintf('processing rhythm: %s-%s\n', tone, rhythm);

        %% get EEG data

        % load data
        datasets = []; 
        
        for i_sub=1:length(subjects)
            
            subject = subjects(i_sub); 
            
            fname = sprintf('AB %s %s P%03d*.lw6', tone, rhythm, subject); 
            
            fprintf('preprocessing: %s\n', fname); 
            
            d = dir(fullfile(load_path, 'derivatives', 'preprocessed_AB', fname)); 
            
            [header, data] =  CLW_load(fullfile(d.folder, d.name)); 

            % low-pass filter 
            [header, data] = RLW_butterworth_filter(header, data, ...
                'filter_type', 'lowpass', 'high_cutoff', 30, 'filter_order', 2); 

            % rereference to mastoid (E57, E100)
            [header, data] = RLW_rereference(header, data,...
                                    'apply_list', {header.chanlocs.labels}, ...
                                    'reference_list', ref_chans); 

            % cut trials into chunks before averaging ? 
            if preproc_chunk_do
                [header_ep, data_ep, msg] = RLW_segmentation_chunk(header, data,...
                    'chunk_onset', preproc_chunk_onset, ...
                    'chunk_duration', preproc_chunk_dur, ...
                    'chunk_interval', preproc_chunk_dur); 
                disp(msg); 
                
                par.trial_dur = preproc_chunk_dur; 
            else
                header_ep = header; 
                data_ep = data; 
            end

            % average epochs 
            [header_ep, data_ep] = RLW_average_epochs(header_ep, data_ep); 

            % select channels 
            % ---------------

            if strcmp(roi,'front')

                % pool channels (spectra averaged across 28 frontocentral channels (14 in each
                % hemisphere))
                chan_sel = front_chan'; 

            elseif strcmp(roi,'best')

                % laod data for selected channel 
                fpath = fullfile(deriv_path,...
                                 sprintf('roi-best/sub-%03d', subject)); 
                fname = sprintf('chan_names.csv'); 
                chan_sel = readtable(fullfile(fpath,fname)); 
                chan_sel = chan_sel.label; 

            end
            
            [header_ep, data_ep] = RLW_arrange_channels(header_ep, data_ep, ...
                                                        chan_sel); 

            
            datasets(i_sub).header = header_ep; 
            datasets(i_sub).data = data_ep; 
            
        end
                
        data = cat(1, datasets.data); 
        
        header = datasets(1).header; 
        header.datasize = size(data);       
              
        n_sub = header.datasize(1); 

        n_chan = header.datasize(2); 

        fs = 1/header.xstep; 

        t = [0 : header.datasize(end)-1] * header.xstep + header.xstart; 

        % make sure we don't have lags longer than half trial duration!
        par.lags_meter_rel = ...
            par.lags_meter_rel(par.lags_meter_rel < par.trial_dur/2); 

        par.lags_meter_unrel = ...
            par.lags_meter_unrel(par.lags_meter_unrel < par.trial_dur/2); 
        

        %% load stimulus

        fname = fullfile(load_path, 'Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat'); 

        coch_output = load(fname); % variables: freq, res_all, rowNames

        fs_coch = coch_output.fs;

        if strcmpi(tone, 'high')
            tone_code = 'H'; 
        else
            tone_code = 'L'; 
        end
        if strcmpi(rhythm,'unsync')
            rhythm_code = 'unsyncopated'; 
        else
            rhythm_code = 'syncopated'; 
        end

        row_idx = ~cellfun(@isempty, regexp(coch_output.cond_names, sprintf('^%s_',tone_code)), 'uni',1) & ...
                  ~cellfun(@isempty, regexp(coch_output.cond_names, sprintf('_%s',rhythm_code)), 'uni',1);

        coch = coch_output.hc(row_idx, :); 
        
        % low-pass filter for nicer acf plot 
        [b,a] = butter(2, 30/(fs_coch/2), 'low'); 
        coch = filtfilt(b, a, coch); 
        
        % cut trials into chunks before averaging ? 
        if preproc_chunk_do
            coch = epoch_chunks(coch, fs_coch, preproc_chunk_dur, ...
                                'start', preproc_chunk_onset); 
            coch = mean(ensure_row(coch), 1); 
        end

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
             ap(i_sub, :, 1, 1, 1, :)] = ...
                                get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'response_f0', par.response_f0, ...
                                       'ap_fit_flims', par.ap_fit_flims, ...
                                       'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                                       'keep_band_around_f0_harmonics', par.ap_band_around_harmonics, ...
                                       'verbose', false);
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

par.roi_name = roi; 

fname = sprintf('exp-infant_apFitMethod-%s_onlyHarm-%s_roi-%s_eegIndividual', ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics),...
                par.roi_name); 

% save table
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save data 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 





