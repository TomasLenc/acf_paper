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

rm_id = [1,4,6,11,14];

par.trial_dur = 50.4;  

par.roi_name = 'frontocentral'; 
par.roi_chans = {'F1', 'Fz', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'Cz', 'C2'}; 

par.ref_name = 'all'; 
% par.ref_chans = 'all'; 


%% allocate table

col_names = {
    'subject', 'rhythm', 'tone', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


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


%% allocate

feat_acf_coch = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'contrast_meter_rel', []); 

feat_acf = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'contrast_meter_rel', []); 

feat_acf_subtracted = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'contrast_meter_rel', []); 

feat_fft_coch = struct('z_meter_rel', []); 

feat_fft = struct('z_meter_rel', []); 

feat_fft_subtracted = struct('z_meter_rel', []); 

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

        % if is unsyncopated, shift the ERP by 3 events!!!
        if strcmpi(rhythm_id, 'unsyncopated')
            offset = 3 * 0.2; 
        else
            offset = 0; 
        end

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
            strfind(coch_output.cond_names, [tone_id, '_standard_', rhythm_id])); 

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
        feat_acf_coch(i_rhythm) = get_acf_features(...
                                    acf_coch, lags_coch, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        % get features for the raw spectra                                    
        feat_fft_coch(i_rhythm) = get_fft_features(mX_coch, freq_coch, ...
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
        ap = nan(size(acf)); 
        parfor i_sub=1:size(data, 1)
            fprintf('sub-%02d\n', i_sub); 
            [acf_subtracted(i_sub, :, 1, 1, 1, :), ~, ...
             ap(i_sub, :, 1, 1, 1, :), ~, ~, ~, ~, ...
             optim_exitflag(i_sub, :)] = ...
                                get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'f0_to_ignore', par.f0_to_ignore, ...
                                       'ap_fit_flims', par.ap_fit_flims, ...
                                       'verbose', false);
        end
        
%         chan_idx = find(strcmp({header.chanlocs.labels}, 'Fz'));
%         [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_exitflag] = ...
%                                     get_acf(data(5, 1, 1, 1, 1, :), ...
%                                             fs, ...
%                                            'ap_fit_method', 'irasa', ...
%                                            'plot_diagnostic', true, ...
%                                            'rm_ap', true, ...
%                                            'f0_to_ignore', 1 / 2.4, ...
%                                            'ap_fit_flims', [0.1, 9]);                                   

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

        feat_acf(i_rhythm) = get_acf_features(acf, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_subtracted(i_rhythm) = get_acf_features(acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 

        % get features for the raw spectra                                    
        tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
        feat_fft(i_rhythm).z_meter_rel = tmp.z_meter_rel; 

        feat_fft(i_rhythm).z_snr = get_z_snr(mX, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 

        % get features for the 1/f-subtracted spectra                                    
        feat_fft_subtracted(i_rhythm) = get_fft_features(mX_subtracted, freq, ...
                                               par.freq_meter_rel, par.freq_meter_unrel);


    %     f = figure('color','white','Position', [944 481 173 179]); 
    %     topoplot(mean(feat_fft.z_snr, 1), ...
    %              header.chanlocs,  ...
    %              'style', 'map', ...
    %              'maplimits', [0, 5], ...
    %              'colormap', parula, ...
    %              'gridscale', 256, ...
    %              'electrodes','off'); 
    %     
    %     f = figure('color','white','Position', [944 481 173 179]); 
    %     topoplot(mean(feat_fft.z_meter_rel, 1), ...
    %              header.chanlocs,  ...
    %              'style', 'map', ...
    %              'maplimits', [-1, 1], ...
    %              'colormap', jet, ...
    %              'gridscale', 256, ...
    %              'electrodes','off');          
    %          
    %     f = figure('color','white','Position', [944 481 173 179]); 
    %     topoplot(mean(feat_acf_subtracted.z_meter_rel, 1), ...
    %              header.chanlocs,  ...
    %              'style', 'map', ...
    %              'maplimits', [-1.5, 1.5], ...
    %              'colormap', jet, ...
    %              'gridscale', 256, ...
    %              'electrodes','off');            


        % add features to table
        rows = [...
            num2cell([1 : n_sub]'), ...
            repmat({rhythm_id}, n_sub, 1), ...
            repmat({tone_id}, n_sub, 1), ...
            num2cell(feat_fft(i_rhythm).z_meter_rel), ...
            num2cell(feat_acf(i_rhythm).z_meter_rel), ...
            num2cell(feat_fft_subtracted(i_rhythm).z_meter_rel), ...
            num2cell(feat_acf_subtracted(i_rhythm).z_meter_rel), ...
            repmat({feat_fft_coch(i_rhythm).z_meter_rel}, n_sub, 1), ...
            repmat({feat_acf_coch(i_rhythm).z_meter_rel}, n_sub, 1), ...
            num2cell(feat_fft(i_rhythm).z_snr) ...
            ];

        tbl = [tbl; rows];
        
    end


end
    
%%

% % colors
% cmap_name = 'Set1'; 
% colors = num2cell(brewermap(n_rhythms, cmap_name), 2); 
% 
% feat = RenameField(feat_acf_subtracted, 'z_meter_rel', 'data');
% feat_orig = RenameField(feat_acf_coch, 'z_meter_rel', 'data');
% 
% % assign labels
% [feat.name] = deal(rhythms{:}); 
% [feat_orig.name] = deal(rhythms{:}); 
% % assign colors
% [feat.color] = deal(colors{:}); 
% [feat_orig.color] = deal(colors{:}); 
% 
% feat(1).data = squeeze(mean(feat(1).data, 2)); 
% feat(2).data = squeeze(mean(feat(2).data, 2)); 
% 
% plot_multiple_cond('plot_legend', false, ...
%                   'zero_line', true, ...
%                   'feat', feat, ...
%                   'feat_orig', feat_orig,...
%                   'ylim_quantile_cutoff', ylim_quantile_cutoff); 
% 
% 
% [H, P, CI, STATS] = ttest(feat(1).data, feat_orig(1).data, 'tail', 'right');
% STATS
% P
% 
% [H, P, CI, STATS] = ttest(feat(2).data, feat_orig(2).data, 'tail', 'right');
% STATS
% P
% 


%% save table

fname = sprintf('exp-lowhigh_apFitMethod-%s_roi-%s_eegIndividual', ...
                par.ap_fit_method, par.roi_name); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




