clear 

addpath('lib')

importLetswave()

rhythms = {'unsyncopated', 'syncopated'}; 

cycle_dur = 2.4; 

% maximum acf lag to consider (for plotting and ap estimation)
max_lag = cycle_dur * 2; 

knee = true; 

subtract_ap = true; 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
lags_meter_rel = [0.8]; 
lags_meter_unrel = [0.6, 1.0]; 

% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 

prec = 1e6; 

col_time = repmat(0, 1, 3); 
col_acf = repmat(0, 1, 3);  
col_ap = [224, 117, 29]/255;  
col_meter_rel = [214, 52, 24]/255; 
col_meter_unrel = [45, 114, 224]/255; 


subject = 'NEGZO'; 
chans = {'H''4', 'H''7', 'R''2', 'H''10'}; 


load_path = '/datadisk/Dropbox/projects/Nancy/data/XPNNC_ClassicRhythmsTAP/derivatives/response-LFP_preproc/'; 

fname = sprintf('%s/sub-%s_response-LFP_preproc.lw6', subject, subject);

[header, data] = CLW_load(fullfile(load_path, fname));

i_rhythm = 1; 

for i_rhythm=1:length(rhythms)

    trig_code = sprintf('%s-listen', rhythms{i_rhythm}); 

    for i_chan=1:length(chans)

        [header_ep, data_ep] = RLW_segmentation(header, data, {trig_code}, ...
                                          'x_start', 0, 'x_duration', 40.8); 

        [header_ep, data_ep] = RLW_arrange_channels(header_ep, data_ep, ...
                                                    chans(i_chan)); 

        [header_ep, data_ep] = RLW_average_epochs(header_ep, data_ep);

        fs = 1/header.xstep; 

        % starting index of acf to keep (note we remove the lag 0)
        min_idx = 2; 
        % last index of acf to keep
        max_idx = round(max_lag * fs); 

        % get acf
        fprintf('calculating acf...\n'); 
        [acf_raw, lags_time] = get_acf(data_ep, fs);

        % take only range of interest
        acf_raw = acf_raw(:, :, :, :, :, min_idx:max_idx); 
        lags_time = lags_time(min_idx:max_idx); 

        % estiamte parameters of the aperiodic component
        ap_par = fit_aperiodic(lags_time, squeeze(acf_raw)', knee); 

        % generate continuous waveform for estimated aperiodic component
        ap = []; 
        ap(1,1,1,1,1,:) = aperiodic(ap_par, lags_time, knee);    


        % subtract ap
        if subtract_ap
            acf = acf_raw - ap; 
        else
            acf = acf_raw; 
        end

        feat = get_acf_features(acf, lags_time, ...
                                lags_meter_rel, lags_meter_unrel, ...
                                lags_meter_unrel_left, lags_meter_unrel_right) ;                       


        % sanity check ERP
        [header_chunk, data_chunk] = RLW_segmentation_chunk(header_ep, data_ep, ...
            'chunk_onset', 0, 'chunk_duration', 2.4, 'chunk_interval', 2.4); 
        erp = squeeze(mean(data_chunk, 1))'; 
        t_erp = [0:length(erp)-1]/fs; 

        % plot
        f = figure('color','white', 'position', [673, 525, 1062, 200], ...
            'name', sprintf('sub-%s %s: %s', subject, trig_code, chans{i_chan})); 
        pnl = panel(f); 

        pnl.pack('h', [20, 40, 40]);
        pnl(1).pack('v', 1);
        pnl(2).pack('v', 1);
        pnl(3).pack('v', 1);

        % plot time-domain 
        pnl(1, 1).select(); 
        plot_time(t_erp, erp, col_time); 
        xlim([0, 2.4])
        axis off

        % plot autocorrelation function
        ax = pnl(2, 1).select(); 
        plot_acf(ax, ...
                 squeeze(acf_raw)', ...
                 lags_time, ...
                 'ap', squeeze(ap)', ...
                 'lags_meter_rel', lags_meter_rel, ...
                 'lags_meter_unrel', lags_meter_unrel, ...
                 'prec', prec, ...
                 'col_acf', col_acf, ...
                 'col_ap', col_ap, ...
                 'col_meter_rel', col_meter_rel, ...
                 'col_meter_unrel', col_meter_unrel ...
                 );            
        yticks([])
        axis off

        % plot corrected autocorrelation function
        ax = pnl(3, 1).select(); 
        plot_acf(ax, ...
                 squeeze(acf)', ...
                 lags_time, ...
                 'lags_meter_rel', lags_meter_rel, ...
                 'lags_meter_unrel', lags_meter_unrel, ...
                 'prec', prec, ...
                 'col_acf', col_acf, ...
                 'col_ap', col_ap, ...
                 'col_meter_rel', col_meter_rel, ...
                 'col_meter_unrel', col_meter_unrel, ...
                 'z_meter_rel', feat.z_meter_rel, ...
                 'ratio_meter_rel', feat.ratio_meter_rel, ...
                 'contrast_meter_rel', feat.contrast_meter_rel ...
                 ); 
             
        yticks([])
        axis off



        pnl(1).xlabel('time (s)')
        pnl(1).ylabel('amplitude')

        pnl(2).xlabel('lag (s)')
        pnl(2).ylabel('autocorrelation')

        pnl(3).xlabel('lag (s)')
        pnl(3).ylabel('autocorrelation')

        pnl.fontsize = 12; 

    end
    
end














