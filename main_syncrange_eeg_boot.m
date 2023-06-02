% function main_syncrange_eeg_boot()
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))


%% parameters

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

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
min_lag = par.min_lag;
max_lag = par.max_lag;

% lags of interest
lags_meter_rel = par.lags_meter_rel;
% meter-unrelated lags 
lags_meter_unrel = par.lags_meter_unrel;

% frequencies of interest
freq_meter_rel = par.freq_meter_rel;
freq_meter_unrel = par.freq_meter_unrel;

% number of neighboring bins for SNR subtraction (FFT only)
noise_bins = par.noise_bins;
noise_bins_snr = [3, 13]; 


%%

n_boot = 500;

n_trials_all = [9:-1:1]; 

rhythms = {'31', '26', '19', '37', '42', '6', '17', '41'}; 

cond_type = 'rhythm'; 

n_rhythms = length(rhythms); 

% colors
cmap_name = 'Set1'; 
colors = num2cell(brewermap(n_rhythms, cmap_name), 2); 

%% allocate table

col_names = {
    'n_trials', 'rhythm', 'boot_sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% RUN

for i_n_trials=1:length(n_trials_all)
    
    n_trials = n_trials_all(i_n_trials); 

    %% allocate
    
    if get_acf_feat_from_x
        
        feat_acf_s = struct(...
            'z_meter_rel', [], 'ratio_meter_rel', [], ...
            'contrast_meter_rel', [], 'mean_meter_rel', []); 
    
        feat_acf = struct(...
            'z_meter_rel', [], 'ratio_meter_rel', [], ...
            'contrast_meter_rel', [], 'mean_meter_rel', []); 
    
        feat_acf_subtracted = struct(...
            'z_meter_rel', [], 'ratio_meter_rel', [], ...
            'contrast_meter_rel', [], 'mean_meter_rel', []); 
    
    else
        
        feat_acf_s = struct(...
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
    
    feat_fft_s = struct('z_meter_rel', []); 
    
    feat_fft = struct('z_meter_rel', []); 
    
    feat_fft_subtracted = struct('z_meter_rel', []); 
    
    feat_ap = struct('offset', [], 'exponent', []); 
    
    %% run
    
    for i_rhythm=1:n_rhythms
        
        rhythm_id = rhythms{i_rhythm};
        
        fprintf('processing rhythm %s\n', rhythm_id);
    
        %% load data
        
        fpath_eeg = fullfile(par.eeg_path, 'syncrange'); 
        
        eeg = load(fullfile(fpath_eeg, ...
            sprintf('exp-syncrange_rhythm-%s_nTrials-%d_eeg.mat', ...
                    rhythm_id, n_trials)));
        
        fs = eeg.fs;
        t = [0 : size(eeg.data, 2) - 1] / fs; 
        
        % take grand average across subjects
        eeg_grand = mean(eeg.data, 1); 
        
        %% load stimulus
        
        fpath_stim = '/DATA1/XPSyncRange/ptb/eeg_7set/to_load'; 
        
        d = dir(fullfile(fpath_stim, sprintf('*rhythm%s_*.wav', rhythm_id)));
        
        [s, fs_s] = audioread(fullfile(fpath_stim, d.name));
        s = s(:, 1)'; 
        env = abs(hilbert(s)); 
        
%         % normalize to eeg
%         s = (s - min(s)) ./ (max(s)-min(s)) ;
%         s = s * (max(eeg_grand) - min(eeg_grand)) + min(eeg_grand); 
        
        t_s = [0 : length(s) - 1] / fs_s; 
        
        % figure
        % plot(t_s, s, 'color', [.8, .8, .8]); 
        % hold on
        % plot(t, eeg_grand, 'linew', 2)
        
        %% process stimulus
        
        % get ACF (withuout aperiodic subtraction)
        [acf_s, lags_s, ~, mX_s, freq_s] = get_acf(...
                                   env, fs_s, ...
                                   'normalize_x', normalize_x, ...
                                   'force_x_positive', force_x_positive, ...
                                   'normalize_acf_to_1', normalize_acf_to_1, ...
                                   'normalize_acf_z', normalize_acf_z ...
                                   );    
                               
        % get ACF features
        feat_acf_s(i_rhythm) = get_acf_features(...
                                    acf_s, lags_s, ...
                                    lags_meter_rel, lags_meter_unrel, ...
                                    'normalize_acf', normalize_acf_vals...
                                    );    
        
        % get features for the raw spectra                                    
        feat_fft_s(i_rhythm) = get_fft_features(mX_s, freq_s, ...
                                         freq_meter_rel, freq_meter_unrel); 
        
        
        %% process EEG
        
        n_sub = size(eeg.data, 1); 
        
        x = nan(n_boot, size(eeg.data, 2));
        
        for i_boot=1:n_boot
            idx = randsample(n_sub, n_sub, true);
            x(i_boot, :) = mean(eeg.data(idx, :), 1); 
        end
    
                
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
        [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_exitflag] = ...
                                    get_acf(x, fs, ...
                                           'rm_ap', true, ...
                                           'f0_to_ignore', 1 / 2.4, ...
                                           'min_freq', 0.1, ...
                                           'max_freq', 9, ...
                                           'get_x_norm', true, ...
                                           'normalize_x', normalize_x, ...
                                           'force_x_positive', force_x_positive, ...
                                           'normalize_acf_to_1', normalize_acf_to_1, ...
                                           'normalize_acf_z', normalize_acf_z ...
                                           );      
        if any(~optim_exitflag)
            warning('ap-fit didnt converge %d/%d reps', sum(~optim_exitflag), n_rhythms); 
        end
        
        feat_ap(i_rhythm).offset = cellfun(@(x) x(1), par_ap);            
        feat_ap(i_rhythm).exponent = cellfun(@(x) x(2), par_ap);            
                                       
        % get features
        % ------------
        
        if get_acf_feat_from_x
    
            feat_acf(i_rhythm) = get_acf_features2(x, fs, ...
                                        lags_meter_rel, lags_meter_unrel);    
    
            feat_acf_subtracted(i_rhythm) = get_acf_features2(x_subtr, fs, ...
                                         lags_meter_rel, lags_meter_unrel); 
            
        else
                    
            feat_acf(i_rhythm) = get_acf_features(acf, lags, ...
                                        lags_meter_rel, lags_meter_unrel, ...
                                        'normalize_acf', normalize_acf_vals);    
    
            feat_acf_subtracted(i_rhythm) = get_acf_features(acf_subtracted, lags, ...
                                         lags_meter_rel, lags_meter_unrel, ...
                                         'normalize_acf', normalize_acf_vals); 
        end
                                     
        % get features for the raw spectra                                    
        tmp = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
        feat_fft(i_rhythm).z_meter_rel = tmp.z_meter_rel; 
                                            
        feat_fft(i_rhythm).z_snr = get_z_snr(mX, freq, par.frex, ...
                                           noise_bins_snr(1), ...
                                           noise_bins_snr(2)); 
    
        % get features for the 1/f-subtracted spectra                                    
        feat_fft_subtracted(i_rhythm) = get_fft_features(mX_subtracted, freq, ...
                                               freq_meter_rel, freq_meter_unrel);
        
                                 
        % plot example 
        % ------------
        
        if plot_example_fig      
            if i_rhythm==1
                f = figure('color','white', ...
                           'position', [95, 67, 1062, 170 * n_rhythms]); 
                pnl_example = panel(f); 
                pnl_example.pack('v', n_rhythms); 
                pnl_example.margin = [5, 10, 25, 25]; 
            end
            rep_to_plot_idx = 1; 
            plot_example(x(rep_to_plot_idx, :), t, ...
                             acf(rep_to_plot_idx, :), lags, ...
                             ap(rep_to_plot_idx, :), ...
                             mX(rep_to_plot_idx, :), freq, ...
                             lags_meter_rel, lags_meter_unrel, ...
                             freq_meter_rel, freq_meter_unrel, ...
                             'pnl', pnl_example(i_rhythm), ...
                             'subplot_proportions', [50, 17, 33], ...
                             'max_lag', max_lag, ...
                             'plot_time_xaxis', i_rhythm == n_rhythms, ...
                             'plot_xlabels', i_rhythm == n_rhythms, ...
                             'plot_xticks', i_rhythm == n_rhythms, ...
                             'plot_features', false, ...
                             'min_lag', 0.2, ...
                             'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                             'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                             'time_col', colors{i_rhythm}, ...
                             'prec', 1e6, ...
                             'fontsize', par.fontsize, ...
                             'normalize_acf_for_plotting', false);                                        
            f.Name = rhythms{i_rhythm};     
            pnl_example(i_rhythm).margintop = 25; 
        end
    
    end
        
    if par.save_figs
       fname = sprintf('exp-syncrange_response-eeg_nTrials-%d_examples.svg', n_trials); 
       save_fig(f, fname)
    end
    
    % assign labels
    [feat_acf.name] = deal(rhythms{:}); 
    [feat_acf_subtracted.name] = deal(rhythms{:}); 
    [feat_fft.name] = deal(rhythms{:}); 
    [feat_fft_subtracted.name] = deal(rhythms{:}); 
    [feat_ap.name] = deal(rhythms{:}); 
    
    % assign colors
    [feat_acf.color] = deal(colors{:}); 
    [feat_fft.color] = deal(colors{:}); 

    [feat_acf_subtracted.color] = deal(colors{:}); 
    [feat_fft_subtracted.color] = deal(colors{:}); 

    [feat_ap.color] = deal(colors{:}); 

    [feat_acf_s.color] = deal([0 0 0]); 
    [feat_fft_s.color] = deal([0 0 0]); 

    %% add features to table
    
    for i_rhythm=1:n_rhythms
    
        rows = [...
            repmat({n_trials}, n_boot, 1), ...
            repmat({feat_acf(i_rhythm).name}, n_boot, 1), ...
            num2cell([1:n_boot]'), ...
            num2cell(feat_fft(i_rhythm).z_meter_rel), ...
            num2cell(feat_acf(i_rhythm).z_meter_rel), ...
            num2cell(feat_fft_subtracted(i_rhythm).z_meter_rel), ...
            num2cell(feat_acf_subtracted(i_rhythm).z_meter_rel), ...
            repmat({feat_fft_s(i_rhythm).z_meter_rel}, n_boot, 1), ...
            repmat({feat_acf_s(i_rhythm).z_meter_rel}, n_boot, 1), ...
            num2cell(feat_fft(i_rhythm).z_snr), ...
            num2cell(feat_ap(i_rhythm).offset), ...
            num2cell(feat_ap(i_rhythm).exponent) ...
            ];
        
        tbl = [tbl; rows];
        
    end
  
    
    
    %% plot
    
    % plot 
    % ----
    if get_acf_feat_from_x
        cond_to_plot = {
            'acf-z_meter_rel'
            'fft-z_meter_rel'
            'fft-z_snr'
            'ap-offset'
            'ap-exponent'
            }; 
    else
        cond_to_plot = {
            'acf-z_meter_rel'
            'fft-z_meter_rel'
            'fft-z_snr'
            'ap-offset'
            'ap-exponent'
            }; 
    end
    
    for i_rhythm=1:length(cond_to_plot)
        
        ytick_at_means = false;
        yaxis_right = false;
        
        switch cond_to_plot{i_rhythm}
            
            case 'acf-mean_meter_rel'
                feat_raw = feat_acf; 
                feat_subtracted = feat_acf_subtracted; 
                feat_sound = feat_acf_s; 
                feat_fieldname = 'mean_meter_rel'; 
                feat_label = 'mean'; 
                tit = 'ACF'; 
            case 'acf-ratio_meter_rel'
                feat_raw = feat_acf; 
                feat_subtracted = feat_acf_subtracted; 
                feat_sound = feat_acf_s; 
                feat_fieldname = 'ratio_meter_rel'; 
                feat_label = 'ratio'; 
                tit = 'ACF'; 
            case 'acf-z_meter_rel'
                feat_raw = feat_acf; 
                feat_subtracted = feat_acf_subtracted; 
                feat_orig = feat_acf_s; 
                feat_fieldname = 'z_meter_rel'; 
                feat_label = 'zscore'; 
                tit = 'ACF'; 
            case 'fft-z_meter_rel'
                feat_raw = feat_fft; 
                feat_subtracted = feat_fft_subtracted; 
                feat_orig = feat_fft_s; 
                feat_fieldname = 'z_meter_rel'; 
                feat_label = 'zscore'; 
                tit = 'FFT'; 
            case 'fft-z_snr'
                feat_raw = feat_fft; 
                feat_subtracted = []; 
                feat_sound = [];
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
        
        f = figure('color', 'white', 'position', [1442 521 350 373]); 
        pnl = panel(f); 
        pnl.pack('h', 2); 
        pnl(1).pack({[0, 0, 1, 1]}); 
        pnl(2).pack({[0, 0, 1, 1]}); 
        
        
        % raw
        ax = pnl(1, 1).select(); 
    
        feat = RenameField(feat_raw, feat_fieldname, 'data');
        feat_orig = RenameField(feat_orig, feat_fieldname, 'data');
    
        plot_multiple_cond('ax', ax, ...
                          'plot_legend', true, ...
                          'feat', feat, ...
                          'feat_orig', feat_orig, ...
                          'ytick_at_means', ytick_at_means, ...
                          'prec', 3, ...
                          'ylim_quantile_cutoff', ylim_quantile_cutoff, ...
                          'point_alpha', 0.2); 
    
        pnl(1).ylabel(sprintf('%s raw', feat_label)); 
        
        if yaxis_right
            ax.YAxisLocation = 'right';
        end
        
        % subtracted
        ax = pnl(2, 1).select(); 
    
        feat = RenameField(feat_subtracted, feat_fieldname, 'data');    
        feat_orig = RenameField(feat_orig, feat_fieldname, 'data');
    
        plot_multiple_cond('ax', ax, ...
                          'plot_legend', false, ...
                          'feat', feat, ...
                          'feat_orig', feat_orig, ...
                          'ytick_at_means', ytick_at_means, ...
                          'prec', 3, ...
                          'ylim_quantile_cutoff', ylim_quantile_cutoff, ...
                          'point_alpha', 0.2); 
    
        pnl(2).ylabel(sprintf('%s 1/f subtr', feat_label)); 
    
        if yaxis_right
            ax.YAxisLocation = 'right';
        end
    
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
        yticks = [];
        c = 0;
        for i_ax=1:length(pnl.children)
            ax = pnl(i_ax, 1).select(); 
            if ~isempty(ax.Children)
                ylims(1) = min(ceil(ax.YLim(1)*prec)/prec, ylims(1)); 
                ylims(2) = max(floor(ax.YLim(2)*prec)/prec, ylims(2)); 
                if isempty(yticks)
                    yticks = ax.YTick;
                else
                    yticks = yticks + ax.YTick;
                end
                c = c+1;
            end
        end
        yticks = yticks ./ c;
        
        for i_ax=1:length(pnl.children)
            ax = pnl(i_ax, 1).select(); 
            if ylims(1) < ylims(2)
                ax.YLim = ylims; 
                ax.YTick = yticks; 
            end
        end
        
        if par.save_figs
           fname = sprintf('exp-syncrange_response-eeg_nTrials-%d_%s_%s.svg', ...
                            n_trials, tit, feat_label);  
           save_fig(f, fname)
        end
           
    end


end

close all


%% save table

fname = sprintf('exp-syncrange_eeg'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 