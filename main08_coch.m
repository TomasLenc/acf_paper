function main08_coch(par)

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))

%% 

tracks = {'H_unsyncopated', 'H_syncopated'}; 

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.0;

% colors
cmap_name = 'Dark2'; 
colors = num2cell(brewermap(2, cmap_name), 2); 

cond_type = 'response'; 


%% 

n_tracks = length(tracks);

for i_track=1:n_tracks

    track = tracks{i_track}; 
    
    % figure
    f = figure('color','white', ...
               'position',  [135 442 1604 502]); 

    pnl = panel(f); 

    pnl.pack('v', [20, 80]); 

    pnl(2).pack('v', n_tracks); 

    example_subplot_proportions = [35, 10, 55];
    
    % allocate
    feat_acf = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'contrast_meter_rel', []); 

    feat_fft = struct('z_meter_rel', []); 

    cond_labels = {}; 

    for i_cond=1:n_tracks
        
        switch i_cond
            case 1
                cond_labels{i_cond} = 'hilbert'; 
                [s, fs] = audioread(fullfile(...
                    par.coch_data_path, ...
                    sprintf(...
                    '%s_200ms_calib75dB_max_short-term_ff_freq1237Hz_onset10_offset50.wav', track...
                    ))); 
                s = s ./ max(abs(s)); 
                x = abs(hilbert(s')); 
                x = decimate(x, 10); 
                fs = fs / 10; 
            case 2
                cond_labels{i_cond} = 'slaney'; 
                coch = load(fullfile(...
                            par.coch_data_path, ...
                            'Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat'...
                            )); 
                x = coch.hc(strcmp(coch.cond_names, track), :); 
                idx = round(0.010 * coch.fs); 
                x(1:idx) = x(idx); 
                x = decimate(x, 10); 
                fs = coch.fs / 10; 
        end

        t = [0:length(x)-1]/fs; 

        % withuout aperiodic subtraction    
        [acf, lags, ~, mX, freq] = get_acf(x, fs);    

        % get features 
        feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

        tmp = get_fft_features(mX, freq, ...
                    par.freq_meter_rel, par.freq_meter_unrel); 
                
        feat_fft(i_cond).z_meter_rel = tmp.z_meter_rel; 

        plot_example(x, t, ...
                     acf, lags, ...
                     [], ...
                     mX, freq, ...
                     par.lags_meter_rel, par.lags_meter_unrel, ...
                     par.freq_meter_rel, par.freq_meter_unrel, ...
                     'pnl', pnl(2, i_cond), ...
                     'subplot_proportions', example_subplot_proportions, ...
                     'min_lag', par.min_lag, ...
                     'max_lag', par.max_lag, ...
                     'max_freq', par.max_freq_plot, ...
                     'plot_time_xaxis', i_cond == n_tracks, ...
                     'plot_xlabels', i_cond == n_tracks, ...
                     'plot_xticks', i_cond == n_tracks, ...
                     'plot_features', false, ...
                     'time_col', colors{i_cond}, ...
                     'prec', 1e8, ...
                     'fontsize', par.fontsize, ...
                     'normalize_acf_for_plotting', false);      
                 
        f.Name = cond_labels{i_cond};     
        pnl(2, i_cond).margintop = 15; 

    end


    % assign labels
    [feat_acf.name] = deal(cond_labels{:}); 
    [feat_fft.name] = deal(cond_labels{:}); 

    % assign colors
    [feat_acf.color] = deal(colors{:}); 
    [feat_fft.color] = deal(colors{:}); 


    % plot 
    % ----
    cond_to_plot = {
        'acf-z_meter_rel'
        'fft-z_meter_rel'
        }; 

    pnl(1).pack('h', 5); 

    for i_cond=1:length(cond_to_plot)

        ytick_at_means = false;
        yaxis_right = false;
        zero_line = false; 
        symmetric_ylims = false; 

        switch cond_to_plot{i_cond}

            case 'acf-z_meter_rel'
                feat_raw = feat_acf; 
                feat_fieldname = 'z_meter_rel'; 
                feat_label = 'zscore'; 
                tit = 'ACF'; 
                zero_line = true; 
                symmetric_ylims = true; 
            case 'fft-z_meter_rel'
                feat_raw = feat_fft; 
                feat_fieldname = 'z_meter_rel'; 
                feat_label = 'zscore'; 
                tit = 'FFT'; 
                zero_line = true; 
                symmetric_ylims = true; 
        end

        pnl(1, i_cond).pack('h', 2); 
        pnl(1, i_cond, 1).pack({[0, 0, 1, 1]}); 
        pnl(1, i_cond, 2).pack({[0, 0, 1, 1]}); 

        % raw
        ax = pnl(1, i_cond, 1, 1).select(); 

        feat = RenameField(feat_raw, feat_fieldname, 'data');

        plot_multiple_cond('ax', ax, ...
                          'plot_legend', true, ...
                          'zero_line', zero_line, ...
                          'ytick_at_means', ytick_at_means, ...
                          'feat', feat, ...
                          'ylim_quantile_cutoff', ylim_quantile_cutoff); 

        pnl(1, i_cond, 1).ylabel(sprintf('%s raw', feat_label)); 

        if yaxis_right
            ax.YAxisLocation = 'right';
        end

        % make the figure nice 
        pnl(1, i_cond).margin = [19, 5, 5, 40]; 

        pnl(1, i_cond).title(tit); 

        pnl(1, i_cond).fontsize = par.fontsize; 

 
    end


    % fix legend position
    sw = 1;
    for i=1:length(f.Children)
        if strcmpi(f.Children(i).Type, 'legend') 
           if sw
               f.Children(i).Title.String = cond_type; 
               f.Children(i).Position(1) = 0.94; 
               f.Children(i).Position(2) = 0.85; 
               sw = 0;
           else
               f.Children(i).Visible = 'off';
           end
        end
    end

    pnl.margin = [15, 10, 25, 15]; 
    pnl(2).margintop = 15;
    
    fname = sprintf('08_coch_track-%s', track); 

    % save parameters 
    save(fullfile(par.fig_path, [fname, '_par.mat']), 'par'); 

    % save figure
    if par.save_figs
        save_fig(f, fullfile(par.fig_path, fname))
    end


end


