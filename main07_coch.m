
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))

coch_data_folder = './data'; 

%% 

track = 'H_syncopated'; 

% make x mean 0 var 1
normalize_x = true; 

% shift x to be positive before calculating acf
force_x_positive = false; 

% scale whole acf between 0 and 1
normalize_acf_to_1 = false; 

% zscore the whole acf 
normalize_acf_z = false; 

% whether to normalize acf values extracted at lags of interest between 0 and 1
normalize_acf_vals = false; 

% plot an example figure for each condition?
plot_example_fig = true; 

% calculate Pearson correlation from requested lags instead of using the
% non-normalized ACF?
get_acf_feat_from_x = false; 

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.0;

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
lags_meter_rel = [0.4, 0.8]; 
lags_meter_unrel = [0.2, 0.6, 1.0]; 

% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 

freq_meter_rel = [1.25 : 1.25 : 5]; 
freq_meter_unrel = setdiff(1/2.4 * [1:12], freq_meter_rel); 

max_freq_plot = 5.5; 

noise_bins = [2, 5]; 

% colors
cmap_name = 'Dark2'; 
colors = num2cell(brewermap(2, cmap_name), 2); 

fontsize = 14; 

save_figs = false; 

%% 

% allocate
if ~get_acf_feat_from_x
    feat_acf = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
        'contrast_meter_rel', []); 
else
    feat_acf = struct(...
        'mean_meter_rel', [], ...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'contrast_meter_rel', []); 
end

feat_fft = struct('z_meter_rel', []); 

n_cond = 2; 
cond_labels = {}; 

for i_cond=1:n_cond
    
    switch i_cond
        case 1
            cond_labels{i_cond} = 'hilbert'; 
            [s, fs] = audioread(fullfile(...
                coch_data_folder, ...
                sprintf(...
                '%s_200ms_calib75dB_max_short-term_ff_freq1237Hz_onset10_offset50.wav', track...
                ))); 
            x = abs(hilbert(s')); 
            x = decimate(x, 10); 
            fs = fs / 10; 
        case 2
            cond_labels{i_cond} = 'slaney'; 
            coch = load(fullfile(...
                        coch_data_folder, ...
                        'Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat'...
                        )); 
            x = coch.hc(strcmp(coch.cond_names, track), :); 
            idx = round(0.010 * coch.fs); 
            x(1:idx) = x(idx); 
            x = decimate(x, 10); 
            fs = coch.fs / 10; 
    end
    
    t = [0:length(x)-1]/fs; 
    
    disp(cond_labels{i_cond}); 

    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, fs, ...
                                       'normalize_x', normalize_x, ...
                                       'force_x_positive', force_x_positive, ...
                                       'normalize_acf_to_1', normalize_acf_to_1, ...
                                       'normalize_acf_z', normalize_acf_z ...
                                       );    

    if ~get_acf_feat_from_x
        feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                    lags_meter_rel, lags_meter_unrel, ...
                                    'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                    'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                    'normalize_acf', normalize_acf_vals);
    else
        feat_acf(i_cond) = get_acf_features2(x, fs, ...
                                    lags_meter_rel, lags_meter_unrel);
    end
    
    feat_fft(i_cond) = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
         
    if plot_example_fig      
        if i_cond==1
            f = figure('color','white', ...
                       'position', [95 67 1062 400]); 
            pnl_example = panel(f); 
            pnl_example.pack('v', n_cond); 
            pnl_example.margin = [5, 10, 25, 25]; 
        end
        rep_to_plot_idx = 1; 
        plot_example(x(rep_to_plot_idx, :), t, ...
                         acf(rep_to_plot_idx, :), lags, ...
                         [], ...
                         mX(rep_to_plot_idx, :), freq, ...
                         lags_meter_rel, lags_meter_unrel, ...
                         freq_meter_rel, freq_meter_unrel, ...
                         'pnl', pnl_example(i_cond), ...
                         'plot_time_xaxis', i_cond == n_cond, ...
                         'plot_xlabels', i_cond == n_cond, ...
                         'plot_xticks', i_cond == n_cond, ...
                         'time_col', colors{i_cond}, ...
                         'prec', 1e6, ...
                         'fontsize', fontsize, ...
                         'normalize_acf_for_plotting', false);                                        
        f.Name = cond_labels{i_cond};  
        pnl_example(i_cond).margintop = 25; 
    end
    
end

if save_figs
   fname = sprintf('figures/07_coch_track-%s_examples.svg', track); 
   print(fname, '-dsvg', '-painters', f);  
end


% assign labels
[feat_acf.name] = deal(cond_labels{:}); 
[feat_fft.name] = deal(cond_labels{:}); 

% assign colors
[feat_acf.color] = deal(colors{:}); 
[feat_fft.color] = deal(colors{:}); 


%%

% plot 
% ----


for i_cond=[2, 4]
    
    switch i_cond
        
        case 1
            feat_raw = feat_acf; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 2
            feat_raw = feat_acf; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'ACF'; 
        case 3
            if ~get_acf_feat_from_x; continue; end
            feat_raw = feat_acf; 
            feat_fieldname = 'mean_meter_rel'; 
            feat_label = 'mean'; 
            tit = 'ACF'; 
        case 4
            feat_raw = feat_fft; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT';     
    end

    f = figure('color', 'white', 'position', [1049 639 232 261]); 
    pnl = panel(f); 
    pnl.pack('h', 2); 
    pnl(1).pack({[0, 0, 1, 1]}); 
    pnl(2).pack({[0, 0, 1, 1]}); 
    
    % raw
    ax = pnl(1, 1).select(); 

    feat = RenameField(feat_raw, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', true, ...
                      'feat', feat, ...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1).ylabel(sprintf('%s raw', feat_label)); 

    % make the figure nice 
    pnl.margin = [19, 5, 5, 40]; 

    pnl.title(tit); 

    pnl.fontsize = fontsize; 

    % fix legend position
    for i=1:length(f.Children)
        if strcmpi(f.Children(i).Type, 'legend')
           f.Children(i).Title.String = 'signal'; 
           f.Children(i).Position(1) = 0; 
           f.Children(i).Position(2) = 0.7; 
        end
    end
    
    ax = pnl(1, 1).select(); 
    prec = 1000; 
    ax.YLim = [-1.3, 1.3];
    ax.YTick = [-1.3, 1.3];
    
    if save_figs
        saveas(f, sprintf('figures/07_coch_track-%s_%s_%s.svg', ...
                            track, tit, feat_label));  
    end
    
end




