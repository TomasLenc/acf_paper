
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))


%% simulate

fs = 200; 

pat = [1 0 1 1 1 1 0 1 1 1 0 0]; % [1 0 1 1 1 1 0 1 1 1 0 0]  [1 1 1 0 1 1 1 0 1 1 0 0]
    
n_cycles = 16; 

grid_ioi = 0.2; 

noise_exponent = -1.5; 

fit_knee = false; 

ir_type = 'square'; 

% number of simulated repetitions 
n_rep = 100; 

% Set true and the signal will be epoched into 1-cycle long chunks and
% averaged. 
do_chunk = true; 
chunk_duration = length(pat) * grid_ioi; 

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
ylim_quantile_cutoff = 0.04; 

% plot an example figure for each condition?
plot_example_fig = true; 

% ------------------------------------------------
cond_type = 'SNR'; 

% snrs = logspace(log10(10), log10(0.2), 6); 
snrs = linspace(1, 0.2, 6); 
% snrs = [0.6, 0.3]; 
% ------------------------------------------------
min_lag = 0; 
max_lag = (grid_ioi * length(pat) * n_cycles) / 2; 


% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
max_lag = 1.21; 
lags_meter_rel = [0.4, 0.8]; 
lags_meter_unrel = [0.2, 0.6, 1.0]; 



% % meter-related lags 
% % ------------------
% 
% % Make sure there's no overlap with muiltiples of meter-unrelated lags, and 
% % also the pattern repetition period. 
% lags_meter_rel = get_lag_harmonics(...
%                             0.8, ...
%                             max_lag,...
%                             'lag_harm_to_exclude', [0.6, 1.0, 2.4] ...
%                             ); 
% 
% % meter-unrelated lags 
% % --------------------
% 
% % Make sure there's no overlap with muiltiples of meter-related lags 
% % (even 0.4 seconds!), and also the pattern repetition period. 
% 
% lags_meter_unrel = get_lag_harmonics(...
%                             0.2, ...
%                             max_lag,...
%                             'lag_harm_to_exclude', [0.4] ...
%                             ); 
% 
% % now take the highest N of each (1/f noise should have little influence at
% % high lags?)
% n_lagz_choose = 2; 
% 
% lags_meter_rel = lags_meter_rel(end-n_lagz_choose+1 : end); 
% 
% [~, idx] = sort(min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel)), [], 1)); 
% 
% lags_meter_unrel = lags_meter_unrel(idx(1:n_lagz_choose)); 
% 
% min_lag = min([lags_meter_rel, lags_meter_unrel]) * 0.9; 
% 
% % make sure one more time that there's no overlap between meter-rel and -unrel !!! 
% assert(~any( min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel))) < 1e-9 ))
% 


% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 

freq_meter_rel = [1.25 : 1.25 : 5]; 
freq_meter_unrel = setdiff(1/2.4 * [1:12], freq_meter_rel); 

frex = sort([freq_meter_rel, freq_meter_unrel]); 
max_freq_plot = 5.5; 

% noise bins for estimating and subtracting the 1/f component
noise_bins = [2, 5]; 

% noies bins for calculating the SNR of the raw spectra
% (harmonic-snippet-zsocre method as used by Rossion)
noise_bins_snr = [3, 13]; 

n_cond = length(snrs); 

% colors
cmap_name = '-OrRd'; 
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(1:n_cond, :); 

fontsize = 14; 


fig_path = 'figures'; 

save_figs = false; 


%% 


if strcmp(ir_type, 'square')
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
elseif strcmp(ir_type, 'erp')
    ir = get_erp_kernel(fs,...
        'amplitudes', 1,...
        't0s', 0, ...
        'taus', 0.050, ...
        'f0s', 7, ...
        'duration', 0.2 ...
        ); 
elseif strcmp(ir_type, 'erp2')
    ir = get_erp_kernel(fs,...
        'amplitudes', [0.4, 0.75],...
        't0s', [0, 0], ...
        'taus', [0.2, 0.050], ...
        'f0s', [1, 7], ...
        'duration', 0.5 ...
        ); 
end

if do_chunk
    % When operating on a chunked->averaged signal, the spectral resolution
    % becomes very small... To prevent crashing, just set the noise bins to 0. 
    noise_bins = [0, 0]; 
    noise_bins_snr = [0, 0]; 
    % Also adapt the maximum lag for plotting. 
    max_lag = grid_ioi * length(pat) / 2; 
end

%% test performance across SNR levels
                             

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

feat_fft = struct('z_meter_rel', [], 'z_snr', []); 

feat_fft_subtracted = struct('z_meter_rel', []); 

cond_labels = {}; 

for i_cond=1:n_cond
    
    snr = snrs(i_cond); 
    cond_labels{i_cond} = sprintf('%.2g', snr); 

    fprintf('calculating snr %d/%d\n', i_cond, n_cond)

    % make clean signal for the whole trial 
    [x_clean, t] = get_s(...
                        pat, ...
                        grid_ioi, ...
                        fs, ...
                        'n_cycles', n_cycles, ...
                        'ir', ir ...
                        );
    
    % generate noisy signal (simulataneously for all repetitions)
    noise = get_colored_noise2([n_rep, length(x_clean)], fs, noise_exponent); 

    % scale the noise to the correct SNR 
    x_clean_rms = rms(x_clean); 
    noise_rms = rms(noise, ndims(noise)); 
    noise_gain = (x_clean_rms ./ noise_rms) / snr; 
    noise = noise .* noise_gain; 
    
    % add signal and noise
    x = x_clean + noise; 
    
    if do_chunk
        
        x_clean_chunked = epoch_chunks(x_clean, fs, chunk_duration); 
        x_chunked = epoch_chunks(x, fs, chunk_duration); 
       
        x_clean = mean(x_clean_chunked, 1); 
        x = squeeze(mean(x_chunked, 1)); 
        
        t = t(1 : length(x_clean)); 
    end
    
    % get acf
    % -------
    
    % clean signal
    [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, fs, ...
                                       'normalize_x', normalize_x, ...
                                       'force_x_positive', force_x_positive, ...
                                       'normalize_acf_to_1', normalize_acf_to_1, ...
                                       'normalize_acf_z', normalize_acf_z ...
                                       );    
                                   
    % withuout aperiodic subtraction    
    [acf, ~, ~, mX, ~] = get_acf(x, fs, ...
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
                                       'get_x_norm', true, ...
                                       'normalize_x', normalize_x, ...
                                       'force_x_positive', force_x_positive, ...
                                       'normalize_acf_to_1', normalize_acf_to_1, ...
                                       'normalize_acf_z', normalize_acf_z ...
                                       );    
                             
                                   
    % get features
    % ------------
    
    if get_acf_feat_from_x
        
        feat_acf_orig(i_cond) = get_acf_features2(x_clean, fs, ...
                                     lags_meter_rel, lags_meter_unrel);         

        feat_acf(i_cond) = get_acf_features2(x, fs, ...
                                    lags_meter_rel, lags_meter_unrel);    

        feat_acf_subtracted(i_cond) = get_acf_features2(x_subtr, fs, ...
                                     lags_meter_rel, lags_meter_unrel); 
        
    else
        
        feat_acf_orig(i_cond) = get_acf_features(acf_clean, lags, ...
                                     lags_meter_rel, lags_meter_unrel, ...
                                     'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                     'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                     'normalize_acf', normalize_acf_vals);         
        
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
            
    % get features for the clean spectra
    feat_fft_orig(i_cond) = get_fft_features(mX_clean, freq,...
                                            freq_meter_rel, freq_meter_unrel); 
    
    % get features for the raw spectra                                    
    tmp = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
    feat_fft(i_cond).z_meter_rel = tmp.z_meter_rel; 
                                        
    feat_fft(i_cond).z_snr = get_z_snr(mX, freq, frex, ...
                                       noise_bins_snr(1), ...
                                       noise_bins_snr(2)); 

    % get features for the 1/f-subtracted spectra                                    
    feat_fft_subtracted(i_cond) = get_fft_features(mX_subtracted, freq, ...
                                           freq_meter_rel, freq_meter_unrel);
    
    % plot example 
    % ------------
    
    if plot_example_fig      
        if i_cond==1
            f = figure('color','white', ...
                       'position', [95 67 1062 170 * n_cond]); 
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
                         'min_lag', min_lag, ...
                         'max_lag', max_lag, ...
                         'plot_time_xaxis', i_cond == n_cond, ...
                         'plot_xlabels', i_cond == n_cond, ...
                         'plot_xticks', i_cond == n_cond, ...
                         'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                         'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                         'time_col', colors{i_cond}, ...
                         'prec', 1e6, ...
                         'fontsize', fontsize, ...
                         'normalize_acf_for_plotting', false);                                        
        f.Name = cond_labels{i_cond};     
        pnl_example(i_cond).margintop = 25; 
    end

end

if save_figs
   fname = sprintf('04_snr_irType-%s_exp-%.1f_nrep-%d_examples.svg', ...
                   ir_type, noise_exponent, n_rep); 
   print(fullfile(fig_path, fname), '-dsvg', '-painters', f);  
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
    cond_to_plot = [1,  3, 4, 5]; 
else
    cond_to_plot = [3, 4, 5]; 
end

for i_cond=cond_to_plot
    
    switch i_cond
        
        case 1
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'mean_meter_rel'; 
            feat_label = 'mean'; 
            tit = 'ACF'; 
        case 2
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 3
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'ACF'; 
        case 4
            feat_raw = feat_fft; 
            feat_subtracted = feat_fft_subtracted; 
            feat_orig = feat_fft_orig; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT'; 
        case 5
            feat_raw = feat_fft; 
            feat_subtracted = []; 
            feat_orig = []; 
            feat_fieldname = 'z_snr'; 
            feat_label = 'zSNR'; 
            tit = 'FFT'; 
    end
    
    f = figure('color', 'white', 'position', [1442 521 350 373]); 
    pnl = panel(f); 
    pnl.pack('h', 2); 
    pnl(1).pack({[0, 0, 1, 1]}); 
    pnl(2).pack({[0, 0, 1, 1]}); 
    
    feat_orig = RenameField(feat_orig, feat_fieldname, 'data');
    
    % raw
    ax = pnl(1, 1).select(); 

    feat = RenameField(feat_raw, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', true, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig,...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1).ylabel(sprintf('%s raw', feat_label)); 

    % subtracted
    ax = pnl(2, 1).select(); 

    feat = RenameField(feat_subtracted, feat_fieldname, 'data');
    
    plot_multiple_cond('ax', ax, ...
                      'plot_legend', false, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig,...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(2).ylabel(sprintf('%s 1/f subtr', feat_label)); 


    % make the figure nice 
    pnl.margin = [19, 5, 5, 40]; 

    pnl.title(tit); 

    pnl.fontsize = fontsize; 

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
        saveas(f, fullfile(fig_path, ...
                           sprintf('04_snr_irType-%s_exp-%.1f_nrep-%d_%s_%s.svg', ...
                                    ir_type, noise_exponent, n_rep, tit, feat_label)));  
    end
   
end

%% t-test actoss first two conditions

if n_cond == 2
    
    [h, p, ci, stats] = ttest(feat_fft_subtracted(1).z_meter_rel, feat_fft_subtracted(2).z_meter_rel); 
    fprintf('\nFFTsubtr z: t(%d) = %.2f, p = %.3f\n', stats.df, stats.tstat, p); 

    [h, p, ci, stats] = ttest(feat_acf_subtracted(1).z_meter_rel, feat_acf_subtracted(2).z_meter_rel); 
    fprintf('\nACFsubtr z: t(%d) = %.2f, p = %.3f\n', stats.df, stats.tstat, p); 
    
end

