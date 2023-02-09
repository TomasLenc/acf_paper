
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))


%% load data

save_figs = true; 


[header, data] = CLW_load(fullfile(...
    par.data_path, ...
    'icfilt(2,8) ep but ds butLP64 sub-005_task-ComplexToneSnr_date-202106171802'...
    )); 
%     'icfilt(2) ep but ds butLP64 sub-004_task-ComplexToneSnr_date-202106171458'...

[header, data] = RLW_rereference(header, data, ...
                        'apply_list', {header.chanlocs.labels}, ...
                        'reference_list',{'TP9', 'TP10'}); 

[header, data] = RLW_arrange_channels(header, data,  ...
                        {'F1','Fz','F2','FC1','FCz','FC2','C1','Cz','C2',}); 

[header, data] = RLW_butterworth_filter(header, data, ...
                                        'filter_type', 'lowpass', ...
                                        'high_cutoff', 20, ...
                                        'filter_order', 2); 
                                    
[header, data] = RLW_segmentation(header, data, {'1'}, 'x_start', 0, 'x_duration', 60); 

fs = 1/header.xstep; 

t = [0 : header.datasize(end)-1]/fs; 


%%

fit_knee = false; 

% make x mean 0 var 1
normalize_x = false; 

% shift x to be positive before calculating acf
force_x_positive = false; 

% scale whole acf between 0 and 1
normalize_acf_to_1 = false; 

% zscore the whole acf 
normalize_acf_z = false; 

% whether to normalize acf values extracted at lags of interest between 0 and 1
normalize_acf_vals = false; 

% calculate Pearson correlation from requested lags instead of using the
% non-normalized ACF?
get_acf_feat_from_x = false; 

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.05; 

% plot an example figure for each condition?
plot_example_fig = true; 

% ------------------------------------------------
cond_type = 'n trials'; 

n_trials = [10, 15, 20, 25, 30]; 

with_replacement = false; 
% ------------------------------------------------

n_cond = length(n_trials); 

% number of repetitions 
n_rep = 50; 


%% 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
min_lag = par.min_lag;
max_lag = par.max_lag;

% meter-related lags 
lags_meter_rel = par.lags_meter_rel;
% meter-unrelated lags 
lags_meter_unrel = par.lags_meter_unrel;

% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = par.lags_meter_unrel_left;
lags_meter_unrel_right = par.lags_meter_unrel_right;

freq_meter_rel = par.freq_meter_rel;
freq_meter_unrel = par.freq_meter_unrel;

noise_bins = par.noise_bins;

%%


% colors
cmap_name = 'PuBu'; 
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(end-n_cond+1:end, :); 



%% 

% allocate
if get_acf_feat_from_x
    
    feat_acf_orig = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'contrast_meter_rel', [], 'mean_meter_rel', []); 

    feat_acf = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'contrast_meter_rel', [], 'mean_meter_rel', []); 

    feat_acf_subtracted = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'contrast_meter_rel', [], 'mean_meter_rel', []); 

else
    
    feat_acf_orig = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
        'contrast_meter_rel', []); 

    feat_acf = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
        'contrast_meter_rel', []); 

    feat_acf_subtracted = struct(...
        'z_meter_rel', [], 'ratio_meter_rel', [], ...
        'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
        'contrast_meter_rel', []); 
end

feat_fft_orig = struct('z_meter_rel', []); 

feat_fft = struct('z_meter_rel', []); 

feat_fft_subtracted = struct('z_meter_rel', []); 

cond_labels = {}; 

for i_cond=1:n_cond
    
    x = nan(n_rep, header.datasize(end)); 
    for i_rep=1:n_rep
        trial_idx = randsample(header.datasize(1), n_trials(i_cond), with_replacement); 
        x(i_rep, :) = squeeze(mean(mean(data(trial_idx, :, 1,1,1, :), 2), 1))'; 
    end
            
    cond_labels{i_cond} = sprintf('%g', n_trials(i_cond)); 

    fprintf('calculating %d/%d\n', i_cond, n_cond)

    % get acf
    % -------
                                       
    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, fs, ...
                               'normalize_x', normalize_x, ...
                               'force_x_positive', force_x_positive, ...
                               'normalize_acf_to_1', normalize_acf_to_1, ...
                               'normalize_acf_z', normalize_acf_z ...
                               );    
                                   
    mX_subtracted = subtract_noise_bins(mX, noise_bins(1),  noise_bins(2)); 
    
    % with aperiodic subtraction    
    [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr] = get_acf(x, fs, ...
                                       'rm_ap', true, ...
                                       'f0_to_ignore', 1/2.4, ...
                                       'min_freq', 0.1, ...
                                       'max_freq', 8, ...
                                       'get_x_norm', true, ...
                                       'normalize_x', normalize_x, ...
                                       'force_x_positive', force_x_positive, ...
                                       'normalize_acf_to_1', normalize_acf_to_1, ...
                                       'normalize_acf_z', normalize_acf_z ...
                                       );                                 
                                   
    % get features
    % ------------
    
    if get_acf_feat_from_x

        feat_acf(i_cond) = get_acf_features2(x, fs, ...
                                    lags_meter_rel, lags_meter_unrel);    

        feat_acf_subtracted(i_cond) = get_acf_features2(x_subtr, fs, ...
                                     lags_meter_rel, lags_meter_unrel); 
        
    else
                
        feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                    lags_meter_rel, lags_meter_unrel, ...
                                    'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                    'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                    'normalize_acf', normalize_acf_vals);    

        feat_acf_subtracted(i_cond) = get_acf_features(acf_subtracted, lags, ...
                                     lags_meter_rel, lags_meter_unrel, ...
                                     'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                     'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                     'normalize_acf', normalize_acf_vals); 
    end
                                 
    feat_fft(i_cond) = get_fft_features(mX, freq, ...
                                            freq_meter_rel, freq_meter_unrel); 

    feat_fft_subtracted(i_cond) = get_fft_features(mX_subtracted, freq, ...
                                           freq_meter_rel, freq_meter_unrel);
    
                             
    % plot example 
    % ------------
    
    if plot_example_fig      
        if i_cond==1
            f = figure('color','white', ...
                       'position', [95, 67, 1062, 150 * n_cond]); 
            pnl_example = panel(f); 
            pnl_example.pack('v', n_cond); 
            pnl_example.margin = [5, 10, 25, 25]; 
        end
        rep_to_plot_idx = 1; 
        plot_example(x(rep_to_plot_idx, :), t, ...
                         acf(rep_to_plot_idx, :), lags, ...
                         ap(rep_to_plot_idx, :), ...
                         mX(rep_to_plot_idx, :), freq, ...
                         lags_meter_rel, lags_meter_unrel, ...
                         freq_meter_rel, freq_meter_unrel, ...
                         'pnl', pnl_example(i_cond), ...
                         'max_lag', max_lag, ...
                         'plot_time_xaxis', i_cond == n_cond, ...
                         'plot_xlabels', i_cond == n_cond, ...
                         'plot_xticks', i_cond == n_cond, ...
                         'min_lag', 0.2, ...
                         'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                         'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                         'time_col', colors{i_cond}, ...
                         'prec', 1e6, ...
                         'fontsize', par.fontsize, ...
                         'normalize_acf_for_plotting', false);                                        
        f.Name = cond_labels{i_cond};     
        pnl_example(i_cond).margintop = 25; 
    end

end

if save_figs
   fname = sprintf('09_eeg_nrep-%d_examples.svg', n_rep); 
   print(fullfile(par.fig_path, fname), '-dsvg', '-painters', f);  
end


% assign labels
[feat_acf.name] = deal(cond_labels{:}); 
[feat_acf_subtracted.name] = deal(cond_labels{:}); 
[feat_fft.name] = deal(cond_labels{:}); 
[feat_fft_subtracted.name] = deal(cond_labels{:}); 

% assign colors
[feat_acf.color] = deal(colors{:}); 
[feat_acf_subtracted.color] = deal(colors{:}); 
[feat_fft.color] = deal(colors{:}); 
[feat_fft_subtracted.color] = deal(colors{:}); 


%%

% plot 
% ----
if get_acf_feat_from_x
    cond_to_plot = [1, 2, 3, 4]; 
else
    cond_to_plot = [2, 3, 4]; 
end

for i_cond=cond_to_plot
    
    switch i_cond
        
        case 1
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_fieldname = 'mean_meter_rel'; 
            feat_label = 'mean'; 
            tit = 'ACF'; 
        case 2
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 3
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'ACF'; 
        case 4
            feat_raw = feat_fft; 
            feat_subtracted = feat_fft_subtracted; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT'; 
    end
    
    f = figure('color', 'white', 'position', [1442 521 350 373]); 
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

    % subtracted
    ax = pnl(2, 1).select(); 

    feat = RenameField(feat_subtracted, feat_fieldname, 'data');
    
    plot_multiple_cond('ax', ax, ...
                      'plot_legend', false, ...
                      'feat', feat, ...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(2).ylabel(sprintf('%s 1/f subtr', feat_label)); 


    % make the figure nice 
    pnl.margin = [19, 5, 5, 40]; 

    pnl.title(tit); 

    pnl.fontsize = par.fontsize; 

    % fix legend position
    for i=1:length(f.Children)
        if strcmpi(f.Children(i).Type, 'legend')
           f.Children(i).Title.String = cond_type; 
           f.Children(i).Position(1) = 0; 
           f.Children(i).Position(2) = 0.7; 
        end
    end

    % get the same ylims across subplots
    prec = 1000; 
    ylims = [Inf, -Inf]; 
    for i_ax=1:length(pnl.children)
        ax = pnl(i_ax, 1).select(); 
        ylims(1) = min(ceil(ax.YLim(1)*prec)/prec, ylims(1)); 
        ylims(2) = max(floor(ax.YLim(2)*prec)/prec, ylims(2)); 
    end
    for i_ax=1:length(pnl.children)
        ax = pnl(i_ax, 1).select(); 
        if ylims(1) < ylims(2)
            ax.YLim = ylims; 
            ax.YTick = ylims; 
        end
    end
    
    if save_figs
       fname = sprintf('09_eeg_nrep-%d_%s_%s.svg', n_rep, tit, feat_label);  
       saveas(f, fullfile(par.fig_path, fname));  
    end
       
end
