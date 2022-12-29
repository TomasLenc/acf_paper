clear 

restoredefaultpath
addpath(genpath('lib'))
addpath(genpath('../acf_tools/src'))

preproc_path = '/datadisk/Dropbox/Tomas_PhD/XP1/XP1_main/matlab/data/EEG/avgRef_timeAvg_mergedParticipants/preprocessed'; 

rhythms = {'unsyncopated', 'syncopated'}; 
tones = {'L', 'H'}; 

cycle_dur = 2.4; 

rmID = [1,4,6,11,14];

knee = true; 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
lags_meter_rel = [0.8]; 
lags_meter_unrel = [0.6, 1.0]; 

prec = 1e6; 

col_time = repmat(0, 1, 3); 
col_acf = repmat(0, 1, 3);  
col_ap = [224, 117, 29]/255;  
col_meter_rel = [214, 52, 24]/255; 
col_meter_unrel = [45, 114, 224]/255; 


chans = {'F1', 'Fz', 'F2', 'FC1', 'FCz', 'FC2'}; 

%%

i_rhythm = 2; 
i_tone = 1; 

for i_rhythm=1:length(rhythms)
    
    for i_tone=1:length(tones)
        
        fname = sprintf('avgRef_timeAvg_mergedParticipants %s_%s(19)', ...
                        tones{i_tone}, rhythms{i_rhythm}); 
        
        [header, data] = CLW_load(fullfile(preproc_path, fname));
        
        [header, data] = RLW_arrange_epochs(header, data, ...
                                setdiff([1:header.datasize(1)], rmID)); 
        
        [header, data] = RLW_segmentation(header, data, {'trig'}, ...
                                          'x_start', 0, 'x_duration', 50.4); 
        
        [header, data] = RLW_rereference(header, data, ...
                                         {header.chanlocs.labels}, ...
                                         {header.chanlocs.labels}); 
        
        [header, data] = RLW_arrange_channels(header, data, chans); 
        
        fs = 1/header.xstep; 
        chan_idx = cellfun(@(ch) find(strcmpi({header.chanlocs.labels}, ch)), ...
                           chans, 'uni', 1); 
        
        % starting index of acf to keep (note we remove the lag 0)
        min_idx = 2; 
        % last index of acf to keep
        max_idx = round(cycle_dur * fs); 
        
        % get acf
        fprintf('calculating acf...\n'); 
        [acf_raw, lags_time] = get_acf(data, fs);

        % take only range of interest
        acf_raw = acf_raw(:, :, :, :, :, min_idx:max_idx); 
        lags_time = lags_time(min_idx:max_idx); 

        % average channels 
        acf_raw = mean(acf_raw, 2); 
        
        
        ap_par = cell(1, size(data, 1)); 
        ap = nan(size(acf_raw)); 
        for i_sub=1:size(data, 1)
            % estiamte parameters of the aperiodic component
            ap_par{i_sub} = fit_aperiodic...
                (lags_time, ...
                squeeze(acf_raw(i_sub,:,1,1,1,:))', ...
                knee); 
            % generate continuous waveform for estimated aperiodic component
            ap(i_sub,1,1,1,1,:) = aperiodic(ap_par{i_sub}, lags_time, knee);    
        end
        
        x = squeeze(mean(mean(acf_raw(:, :, 1, 1, 1, :), 1), 2)); 
        
        f = figure('color','white', 'position', [668 62 257 900], ...
                   'name', sprintf('%s-%s', rhythms{i_rhythm}, tones{i_tone})); 
        pnl = panel(f); 
        pnl.pack('v', 14); 
        for sub2plot=1:14
            ax = pnl(sub2plot).select(); 
            plot_acf(ax, ...
                     squeeze(acf_raw(sub2plot,1,1,1,1,:))', ...
                     lags_time, ...
                     'ap', squeeze(ap(sub2plot,1,1,1,1,:))', ...
                     'lags_meter_rel', lags_meter_rel, ...
                     'lags_meter_unrel', lags_meter_unrel, ...
                     'prec', prec, ...
                     'col_acf', col_acf, ...
                     'col_ap', col_ap, ...
                     'col_meter_rel', col_meter_rel, ...
                     'col_meter_unrel', col_meter_unrel ...
                     );                  
        end
        
        pnl.de.margin = [1,1,1,1]; 
        pnl.margin = [15, 6, 10, 5]; 
        
        % sanity check ERP
        figure; 
        [header_chunk, data_chunk] = RLW_segmentation_chunk(header, data, ...
            'chunk_onset', 0, 'chunk_duration', 2.4, 'chunk_interval', 2.4); 
        x = squeeze(mean(mean(data_chunk(:, :, 1, 1, 1, :), 1), 2)); 
        t = [0:length(x)-1]/fs; 
        plot(t, x, 'linew', 2)
        
    end
    
end














