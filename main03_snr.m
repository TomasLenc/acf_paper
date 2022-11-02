
clear 
fprintf('\n');
addpath(genpath('lib'))
addpath('/datadisk/Dropbox/projects/autocorrelation/acf_tools/src'); 

%% simulate

fs = 200; 

exponent = -1; 

fit_knee = false; 

normalize_acf = false; 

% cutoff 5% of extreme points for plotting 
ylim_quantile_cutoff = 0.05; 

snrs = logspace(log10(10), log10(0.2), 6); 

% number of simulated repetitions 
n_rep = 100; 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
lags_meter_rel = [0.8]; 
lags_meter_unrel = [0.6, 1.0]; 

% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 


freq_meter_rel = [1.25 : 1.25 : 5]; 
freq_meter_unrel = setdiff(1/2.4 * [1:12], freq_meter_rel); 

max_freq_plot = 5.5; 
noise_bins = [2, 5]; 

% colors
cmap_name = '-OrRd'; 
colors = num2cell(brewermap(length(snrs) + length(snrs), cmap_name), 2); 
colors = colors(1:length(snrs), :); 

%% 

ir_par = []; 
ir_par.type = 'square'; 
ir_par.eventDur = 0.050; % 4 ms for click, 100ms for tone
ir_par.rampon = 0.0;
ir_par.rampoff = 0.0; 
ir = get_ir(ir_par, fs); 

% make whole signal 
[x_clean, t] = get_s(...
                    [1 0 1 1 1 1 0 1 1 1 0 0], ...
                    16, ...
                    0.2, ...
                    0, ...
                    ir, ...
                    fs, ...
                    'emph_magn', 0, ...
                    'emph_period', 4 ...
                    );
                
x_clean = (x_clean - mean(x_clean)) ./ std(x_clean); 

% make the signal strictly positive
x_clean = x_clean - min(x_clean); 

% no noise (clean signal) 
[acf, lags, ap, mX, freq] = get_acf(x_clean, fs); 

feat_acf_orig = get_acf_features(acf, lags, ...
                             lags_meter_rel, lags_meter_unrel, ...
                             'lags_meter_unrel_left', lags_meter_unrel_left, ...
                             'lags_meter_unrel_right', lags_meter_unrel_right, ...
                             'normalize_acf', normalize_acf); 

feat_fft_orig = get_fft_features(mX, freq, ...
                                 freq_meter_rel, freq_meter_unrel);

%% test performance across SNR levels
                             
% allocate
feat_acf = feat_acf_orig; 
feat_acf_subtracted = feat_acf_orig; 
feat_fft = feat_fft_orig; 
feat_fft_subtracted = feat_fft_orig; 
           
cond_labels = {}; 

for i_cond=1:length(snrs)
    
    snr = snrs(i_cond); 
    cond_labels{i_cond} = sprintf('%.2g', snr); 

    fprintf('calculating snr %d/%d\n', i_cond, length(snrs))

    % generate noisy signal (simulataneously for all repetitions)
    noise = get_colored_noise2([n_rep, length(x_clean)], fs, exponent); 

    x = x_clean + noise./snr; 

    % make the signal purely positive
    x = x - min(x, [], ndims(x)); 
    
    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, fs); 

    feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                lags_meter_rel, lags_meter_unrel, ...
                                'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                'normalize_acf', normalize_acf);
    
    feat_fft(i_cond) = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 

    
    % with aperiodic subtraction    
    [acf_subtracted, ~, ap] = get_acf(x, fs, 'rm_ap', true,  'f0_to_ignore', 1/2.4); 

    feat_acf_subtracted(i_cond) = get_acf_features(acf_subtracted, lags, ...
                                     lags_meter_rel, lags_meter_unrel, ...
                                     'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                     'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                     'normalize_acf', normalize_acf); 
    
    mX_subtracted = subtract_noise_bins(mX, noise_bins(1),  noise_bins(2)); 

    feat_fft_subtracted(i_cond) = get_fft_features(mX_subtracted, freq, ...
                                           freq_meter_rel, freq_meter_unrel);

    rep_to_plot_idx = 1; 
    f = plot_example(x(rep_to_plot_idx, :), t, ...
                     acf(rep_to_plot_idx, :), lags, ...
                     ap(rep_to_plot_idx, :), ...
                     mX(rep_to_plot_idx, :), freq, ...
                     lags_meter_rel, lags_meter_unrel, ...
                     freq_meter_rel, freq_meter_unrel, ...
                     'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                     'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                     'time_col', colors{i_cond});                                        
    f.Name = cond_labels{i_cond};     
                                    
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


for i_cond=[1, 3]
    
    switch i_cond
        
        case 1
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 2
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'contrast_meter_rel'; 
            feat_label = 'contrast'; 
            tit = 'ACF'; 
        case 3
            feat_raw = feat_fft; 
            feat_subtracted = feat_fft_subtracted; 
            feat_orig = feat_fft_orig; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT'; 
    end
            

    f = figure('color', 'white', 'position', [945 371 350 373]); 
    pnl = panel(f); 
    pnl.pack('h', 2); 

    % raw
    ax = pnl(1).select(); 

    feat = RenameField(feat_raw, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', true, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig.(feat_fieldname),...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1).ylabel(sprintf('%s raw', feat_label)); 

    % subtracted
    ax = pnl(2).select(); 

    feat = RenameField(feat_subtracted, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', false, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig.(feat_fieldname), ...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(2).ylabel(sprintf('%s 1/f subtr', feat_label)); 


    % make the figure nice 
    pnl.margin = [19, 5, 5, 40]; 

    pnl.title(tit); 

    pnl.fontsize = 12; 

    % fix legend position
    for i=1:length(f.Children)
        if strcmpi(f.Children(i).Type, 'legend')
           f.Children(i).Title.String = 'SNR'; 
           f.Children(i).Position(1) = 0; 
           f.Children(i).Position(2) = 0.7; 
        end
    end

    % get the same ylims across subplots
    prec = 1000; 
    ylims = [Inf, -Inf]; 
    for i_ax=1:length(pnl.descendants)
        ax = pnl(i_ax).select(); 
        ylims(1) = min(ceil(ax.YLim(1)*prec)/prec, ylims(1)); 
        ylims(2) = max(floor(ax.YLim(2)*prec)/prec, ylims(2)); 
    end
    for i_ax=1:length(pnl.descendants)
        ax = pnl(i_ax).select(); 
        ax.YLim = ylims; 
        ax.YTick = ylims; 
    end

end

