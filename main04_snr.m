function main04_snr()
% clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%% simulate

noise_exponent = -1.5; 

fit_knee = false; 

noise_type = 'eeg'; % eeg, fractal

ir_type = 'square'; 

% number of simulated repetitions 
n_rep = 200; 

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.04; 

% ------------------------------------------------
cond_type = 'SNR'; 

if strcmp(noise_type, 'eeg')
    snrs = logspace(log10(0.2), log10(2), 5); 
%     snrs(1) = 1/1e7;
else
    snrs = logspace(log10(0.2), log10(2), 5); 
end
% snrs = linspace(1, 0.2, 6); 
% snrs = [0.6, 0.3]; 
% ------------------------------------------------

%%

n_cond = length(snrs); 

% colors
cmap_name = 'OrRd'; 
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(end-n_cond+1:end, :); 


%% 

if strcmp(ir_type, 'square')
    ir = get_square_kernel(par.fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
elseif strcmp(ir_type, 'erp')
    ir = get_erp_kernel(par.fs,...
        'amplitudes', 1,...
        't0s', 0, ...
        'taus', 0.050, ...
        'f0s', 7, ...
        'duration', 0.2 ...
        ); 
elseif strcmp(ir_type, 'erp2')
    ir = get_erp_kernel(par.fs,...
        'amplitudes', [0.4, 0.75],...
        't0s', [0, 0], ...
        'taus', [0.2, 0.050], ...
        'f0s', [1, 7], ...
        'duration', 0.5 ...
        ); 
end

%% test performance across SNR levels
                            
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

feat_fft_orig = struct('z_meter_rel', []); 

feat_fft = struct('z_meter_rel', [], 'z_snr', []); 

feat_fft_subtracted = struct('z_meter_rel', []); 

feat_ap = struct('offset', [], 'exponent', []); 

cond_labels = {}; 

ymax_mX = -Inf;
ymax_mX_subtracted = -Inf;



%% generate signal

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );

%% genearet noise

if strcmp(noise_type, 'fractal')

    noise = get_colored_noise2([n_rep, length(x_clean)], par.fs, noise_exponent); 

elseif strcmp(noise_type, 'eeg')

    trial_dur = par.n_cycles * length(par.pat) * par.grid_ioi; 

    noise = prepare_eeg_noise(n_rep, trial_dur);    

else
    error('noise type "%s" not implemented', noise_type);
end


%% run

f = figure('color','white', ...
           'position', [163 1222 1604 1127]); 
       
pnl = panel(f); 

pnl.pack('v', [20, 80]); 

pnl(2).pack('v', n_cond); 

example_subplot_proportions = [40, 15, 45];


for i_cond=1:n_cond
    
    snr = snrs(i_cond); 
    cond_labels{i_cond} = sprintf('%.2g', snr); 

    fprintf('calculating snr %d/%d\n', i_cond, n_cond)
    
    % scale the noise to the correct SNR 
    x = add_signal_noise(repmat(x_clean, n_rep, 1), noise, snr);
    
    % get acf
    % -------
    
    % clean signal
    [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    
                                   
    % withuout aperiodic subtraction    
    [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    
                                   
    mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1), par.noise_bins(2)); 
    
    % with aperiodic subtraction    
    [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_flag] = ...
                                get_acf(x, par.fs, ...
                                       'rm_ap', true, ...
                                       'f0_to_ignore', 1/2.4, ...
                                       'min_freq', 0.1, ...
                                       'max_freq', 9);  
                                   
    feat_ap(i_cond).offset = cellfun(@(x) x(1), par_ap);            
    feat_ap(i_cond).exponent = cellfun(@(x) x(2), par_ap);            
    
    % get features
    % ------------
 
    feat_acf_orig(i_cond) = get_acf_features(acf_clean, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

    feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    

    feat_acf_subtracted(i_cond) = get_acf_features(acf_subtracted, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel); 

    % get features for the clean spectra
    feat_fft_orig(i_cond) = get_fft_features(mX_clean, freq,...
                                            par.freq_meter_rel, par.freq_meter_unrel); 
    
    % get features for the raw spectra                                    
    tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
    feat_fft(i_cond).z_meter_rel = tmp.z_meter_rel; 
                                        
    feat_fft(i_cond).z_snr = get_z_snr(mX, freq, par.frex, ...
                                       par.noise_bins_snr(1), ...
                                       par.noise_bins_snr(2)); 

    % get features for the 1/f-subtracted spectra                                    
    feat_fft_subtracted(i_cond) = get_fft_features(mX_subtracted, freq, ...
                                           par.freq_meter_rel, par.freq_meter_unrel);
    
    % plot example 
    % ------------

    rep_to_plot_idx = 1; 

    % update yaxis maximum for FFT
    frex_idx = dsearchn(freq', par.frex');
    amps = mX(rep_to_plot_idx, frex_idx);
    ymax_mX = max(ymax_mX, max(amps));
    amps = mX_subtracted(rep_to_plot_idx, frex_idx);
    ymax_mX_subtracted = max(ymax_mX_subtracted, max(amps));

    plot_example(x(rep_to_plot_idx, :), t, ...
                     acf(rep_to_plot_idx, :), lags, ...
                     ap(rep_to_plot_idx, :), ...
                     mX(rep_to_plot_idx, :), freq, ...
                     par.lags_meter_rel, par.lags_meter_unrel, ...
                     par.freq_meter_rel, par.freq_meter_unrel, ...
                     'pnl', pnl(2, i_cond), ...
                     'subplot_proportions', example_subplot_proportions, ...
                     'min_lag', par.min_lag, ...
                     'max_lag', par.max_lag, ...
                     'max_freq', par.max_freq_plot, ...
                     'plot_time_xaxis', i_cond == n_cond, ...
                     'plot_xlabels', i_cond == n_cond, ...
                     'plot_xticks', i_cond == n_cond, ...
                     'plot_features', false, ...
                     'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                     'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                     'time_col', colors{i_cond}, ...
                     'prec', 1e8, ...
                     'fontsize', par.fontsize, ...
                     'normalize_acf_for_plotting', false);                                        
    f.Name = cond_labels{i_cond};     
    pnl(2, i_cond).margintop = 24; 


end

for i_cond=1:n_cond
    ax = pnl(2, i_cond, 2, 1).select();
    ax.YLim = [0, ymax_mX];
    ax = pnl(2, i_cond, 2, 2).select();
    ax.YLim = [0, ymax_mX_subtracted];
end

% assign labels
[feat_acf.name] = deal(cond_labels{:}); 
[feat_acf_subtracted.name] = deal(cond_labels{:}); 
[feat_fft.name] = deal(cond_labels{:}); 
[feat_fft_subtracted.name] = deal(cond_labels{:}); 
[feat_ap.name] = deal(cond_labels{:}); 

% assign colors
[feat_acf.color] = deal(colors{:}); 
[feat_acf_subtracted.color] = deal(colors{:}); 
[feat_fft.color] = deal(colors{:}); 
[feat_fft_subtracted.color] = deal(colors{:}); 
[feat_ap.color] = deal(colors{:}); 


%%

% plot 
% ----
cond_to_plot = {
    'acf-z_meter_rel'
    'fft-z_meter_rel'
    'fft-z_snr'
    'ap-offset'
    'ap-exponent'
    }; 

pnl(1).pack('h', length(cond_to_plot)); 


for i_cond=1:length(cond_to_plot)
    
    ytick_at_means = false;
    yaxis_right = false;
    zero_line = false; 
    symmetric_ylims = false; 
    
    switch cond_to_plot{i_cond}
        
        case 'acf-mean_meter_rel'
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'mean_meter_rel'; 
            feat_label = 'mean'; 
            tit = 'ACF'; 
        case 'acf-ratio_meter_rel'
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 'acf-z_meter_rel'
            feat_raw = feat_acf; 
            feat_subtracted = feat_acf_subtracted; 
            feat_orig = feat_acf_orig; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'ACF'; 
            zero_line = true; 
            symmetric_ylims = true; 
        case 'fft-z_meter_rel'
            feat_raw = feat_fft; 
            feat_subtracted = feat_fft_subtracted; 
            feat_orig = feat_fft_orig; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT'; 
            zero_line = true; 
            symmetric_ylims = true; 
        case 'fft-z_snr'
            feat_raw = feat_fft; 
            feat_subtracted = []; 
            feat_orig = []; 
            feat_fieldname = 'z_snr'; 
            feat_label = 'zSNR'; 
            tit = 'FFT'; 
            ytick_at_means = true;
            yaxis_right = true;
        case 'ap-offset'
            feat_raw = feat_ap; 
            feat_subtracted = []; 
            feat_orig = []; 
            feat_fieldname = 'offset'; 
            feat_label = 'offset'; 
            tit = 'AP'; 
        case 'ap-exponent'
            feat_raw = feat_ap; 
            feat_subtracted = []; 
            feat_orig = []; 
            feat_fieldname = 'exponent'; 
            feat_label = 'exponent'; 
            tit = 'AP'; 
    end
    
    pnl(1, i_cond).pack('h', 2); 
    pnl(1, i_cond, 1).pack({[0, 0, 1, 1]}); 
    pnl(1, i_cond, 2).pack({[0, 0, 1, 1]}); 
    
    feat_orig = RenameField(feat_orig, feat_fieldname, 'data');
    
    % raw
    ax = pnl(1, i_cond, 1, 1).select(); 

    feat = RenameField(feat_raw, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', true, ...
                      'zero_line', zero_line, ...
                      'ytick_at_means', ytick_at_means, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig,...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1, i_cond, 1).ylabel(sprintf('%s raw', feat_label)); 

    if yaxis_right
        ax.YAxisLocation = 'right';
    end

    % subtracted
    ax = pnl(1, i_cond, 2, 1).select(); 

    feat = RenameField(feat_subtracted, feat_fieldname, 'data');
    
    plot_multiple_cond('ax', ax, ...
                      'plot_legend', false, ...
                      'zero_line', zero_line, ...
                      'ytick_at_means', ytick_at_means, ...
                      'feat', feat, ...
                      'feat_orig', feat_orig,...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1, i_cond, 2).ylabel(sprintf('%s 1/f subtr', feat_label)); 

    if yaxis_right
        ax.YAxisLocation = 'right';
    end
    
    % make the figure nice 
    pnl(1, i_cond).margin = [19, 5, 5, 40]; 

    pnl(1, i_cond).title(tit); 

    pnl(1, i_cond).fontsize = par.fontsize; 

    % get the same ylims across subplots
    prec = 1000; 
    ylims = [Inf, -Inf]; 
    yticks = [];
    c = 0;
    for i_ax=1:length(pnl(1, i_cond).children)
        ax = pnl(1, i_cond, i_ax, 1).select(); 
        if ~isempty(ax.Children)
            ylims(1) = min(ceil(ax.YLim(1)*prec)/prec, ylims(1)); 
            ylims(2) = max(floor(ax.YLim(2)*prec)/prec, ylims(2)); 
            if isempty(yticks)
                yticks = ax.YTick;
            else
                yticks = yticks + ax.YTick;
            end
            c = c+1;
        else
            ax.Visible = 'off';
            pnl(1, i_cond, i_ax).ylabel('');
            pnl(1, i_cond, i_ax).xlabel('');
        end
    end
    yticks = yticks ./ c;
    
    if symmetric_ylims
        ylims = [-max(ylims), +max(ylims)]; 
        yticks = ylims; 
    end
    
    for i_ax=1:length(pnl(1, i_cond).children)
        ax = pnl(1, i_cond, i_ax, 1).select(); 
        if ylims(1) < ylims(2)
            ax.YLim = ylims; 
            ax.YTick = yticks; 
        end
    end
    

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
pnl(2).margintop = 35;

if strcmp(noise_type, 'fractal')
    fname = sprintf('04_snr_irType-%s_exp-%.1f_nrep-%d', ...
                   ir_type, noise_exponent, n_rep); 
elseif strcmp(noise_type, 'eeg')
    fname = sprintf('04_snr_irType-%s_noise-eeg_nrep-%d', ...
                   ir_type, n_rep); 
else
    error('noise type "%s" not implemented', noise_type);
end

if par.save_figs
   save_fig(f, fname)
end

% save parameters 
save(fullfile(par.fig_path, [fname, '_par.mat']), 'par'); 
