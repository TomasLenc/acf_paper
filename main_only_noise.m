% function main04_snr()

addpath(genpath(par.lw_path)); 

par = get_par(); 

%% simulate

noise_exponent = -1.5; 

fit_knee = false; 

noise_type = 'eeg'; % eeg, fractal

% number of simulated repetitions 
n_rep = 1000; 

% Set true and the signal will be epoched into 1-cycle long chunks and
% averaged. 
do_chunk = false; 
chunk_duration = length(par.pat) * par.grid_ioi; 

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
noise_bins_snr = par.noise_bins_snr;


trial_dur = par.n_cycles * length(par.pat) * par.grid_ioi; 
N = round(trial_dur * par.fs); 

if do_chunk
    % When operating on a chunked->averaged signal, the spectral resolution
    % becomes very small... To prevent crashing, just set the noise bins to 0. 
    noise_bins = [0, 0]; 
    noise_bins_snr = [0, 0]; 
    % Also adapt the maximum lag for plotting. 
    max_lag = par.grid_ioi * length(par.pat) / 2; 
end


%%

% generate noisy signal (simulataneously for all repetitions) 
if strcmp(noise_type, 'fractal')

    x = get_colored_noise2([n_rep, N], par.fs, noise_exponent); 

elseif strcmp(noise_type, 'eeg')

    x = prepare_eeg_noise(n_rep, trial_dur); 

else
    error('noise type "%s" not implemented', noise_type);
end


%%

% get acf
% -------

% withuout aperiodic subtraction    
[acf, lags, ~, mX, freq] = get_acf(x, par.fs, ...
                           'normalize_x', normalize_x, ...
                           'force_x_positive', force_x_positive, ...
                           'normalize_acf_to_1', normalize_acf_to_1, ...
                           'normalize_acf_z', normalize_acf_z ...
                           );    

mX_subtracted = subtract_noise_bins(mX, noise_bins(1),  noise_bins(2)); 

% with aperiodic subtraction    
[acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr] = get_acf(x, par.fs, ...
                                   'rm_ap', true, ...
                                   'f0_to_ignore', 1/2.4, ...
                                   'min_freq', 0.1, ...
                                   'max_freq', 9, ...
                                   'get_x_norm', true, ...
                                   'normalize_x', normalize_x, ...
                                   'force_x_positive', force_x_positive, ...
                                   'normalize_acf_to_1', normalize_acf_to_1, ...
                                   'normalize_acf_z', normalize_acf_z ...
                                   );  

feat_ap.offset = cellfun(@(x) x(1), par_ap);            
feat_ap.exponent = cellfun(@(x) x(2), par_ap);            

% get features
% ------------

if get_acf_feat_from_x

    feat_acf = get_acf_features2(x, par.fs, ...
                                lags_meter_rel, lags_meter_unrel);    

    feat_acf_subtracted = get_acf_features2(x_subtr, par.fs, ...
                                 lags_meter_rel, lags_meter_unrel); 

else

    feat_acf = get_acf_features(acf, lags, ...
                                lags_meter_rel, lags_meter_unrel, ...
                                'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                'normalize_acf', normalize_acf_vals);    

    feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                 lags_meter_rel, lags_meter_unrel, ...
                                 'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                 'lags_meter_unrel_right', lags_meter_unrel_right, ...
                                 'normalize_acf', normalize_acf_vals); 
end

% get features for the raw spectra                                    
tmp = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
feat_fft.z_meter_rel = tmp.z_meter_rel; 

feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                   noise_bins_snr(1), ...
                                   noise_bins_snr(2)); 

% get features for the 1/f-subtracted spectra                                    
feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                       freq_meter_rel, freq_meter_unrel);




col = [0, 0, 0];

% assign labels
feat_acf.name = ''; 
feat_acf_subtracted.name = ''; 
feat_fft.name = ''; 
feat_fft_subtracted.name = ''; 
feat_ap.name = ''; 

% assign colors
feat_acf.color = col; 
feat_acf_subtracted.color = col; 
feat_fft.color = col; 
feat_fft_subtracted.color = col; 
feat_ap.color = col; 

%%

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


%%


close all

f = figure('color', 'white'); 
pnl = panel(f); 
pnl.pack('h', length(cond_to_plot)); 



i_subplot = 1;
ax = pnl(i_subplot).select(); 
plot_points(ax, 0, feat_fft.z_meter_rel, ...
            'opacity', 0.4, ...
            'col', feat_fft.color); 
hold(ax, 'on');
plot(ax, ...
     [-0.4, +0.4], ...
     [0, 0], ...
     '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)

ax.XLim = [-0.4, 0.4]; 
ax.XAxis.Visible = 'off';
ax.YTick = [0];
ax.YLim = [-max(abs(ax.YLim)), max(abs(ax.YLim))];
pnl(i_subplot).ylabel('z beat rel');
pnl(i_subplot).title('fft raw'); 


i_subplot = 2;
ax = pnl(i_subplot).select(); 
plot_points(ax, 0, feat_fft_subtracted.z_meter_rel, ...
            'opacity', 0.4, ...
            'col', feat_fft_subtracted.color); 
hold(ax, 'on');
plot(ax, ...
     [-0.4, +0.4], ...
     [0, 0], ...
     '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)
ax.XLim = [-0.4, 0.4]; 
ax.XAxis.Visible = 'off';
ax.YTick = [0];
ax.YLim = [-max(abs(ax.YLim)), max(abs(ax.YLim))];
pnl(i_subplot).ylabel('z beat rel');
pnl(i_subplot).title('fft subtr'); 


i_subplot = 3;
ax = pnl(i_subplot).select(); 
plot_points(ax, 0, feat_acf.z_meter_rel, ...
            'opacity', 0.4, ...
            'col', feat_acf.color); 
hold(ax, 'on');
plot(ax, ...
     [-0.4, +0.4], ...
     [0, 0], ...
     '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)
ax.XLim = [-0.4, 0.4]; 
ax.XAxis.Visible = 'off';
ax.YTick = [0];
ax.YLim = [-max(abs(ax.YLim)), max(abs(ax.YLim))];
pnl(i_subplot).ylabel('z beat rel');
pnl(i_subplot).title('acf raw'); 


i_subplot = 4;
ax = pnl(i_subplot).select(); 
plot_points(ax, 0, feat_acf_subtracted.z_meter_rel, ...
            'opacity', 0.4, ...
            'col', feat_acf.color); 
hold(ax, 'on');
plot(ax, ...
     [-0.4, +0.4], ...
     [0, 0], ...
     '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)

ax.XLim = [-0.4, 0.4]; 
ax.XAxis.Visible = 'off';
ax.YTick = [0];
ax.YLim = [-max(abs(ax.YLim)), max(abs(ax.YLim))];
pnl(i_subplot).ylabel('z beat rel');
pnl(i_subplot).title('acf subtr'); 


pnl.de.margin = [5, 0, 0, 0]; 
pnl.margin = [15, 5, 0, 15]; 

if strcmp(noise_type, 'fractal')
    fname = sprintf('onlyNoise_noise-fractal_exp-%.1f_nrep-%d', ...
                   noise_exponent, n_rep); 
elseif strcmp(noise_type, 'eeg')
    fname = sprintf('onlyNoise_noise-eeg_nrep-%d', ...
                   n_rep); 
else
    error('noise type "%s" not implemented', noise_type);
end

if par.save_figs
   save_fig(f, fname)
end

   

%% Bayes factor t-test against 0 

tbl = [
    num2cell(feat_fft.z_meter_rel), ...
    num2cell(feat_fft_subtracted.z_meter_rel), ...
    num2cell(feat_acf.z_meter_rel), ...
    num2cell(feat_acf_subtracted.z_meter_rel)
    ];

tbl = cell2table(tbl, ...
    'VariableNames', {'fft_raw', 'fft_subtr', 'acf_raw', 'acf_subtr'}); 

writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

%% 

mX_grand = mean(mX, 1); 
% mX_err = std(mX, [], 1) / sqrt(n_rep) * norminv(1 - 0.025); 
mX_err = std(mX, [], 1); 

col = [165, 52, 217] / 255; 

f = figure('color', 'white', 'Position', [769 470 137 132]); 
pnl = panel(f); 

ax = pnl.select(); 
fill(ax, [freq, flip(freq)], [mX_grand+mX_err, flip(mX_grand-mX_err)],...
    col, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
hold(ax, 'on'); 
plot(ax, freq(2:end), mX_grand(2:end), 'linew', 0.5, 'color', col); 
ax.XLim = [0.1, 15]; 
ax.XTick = [ax.XLim(1), ax.XLim(end)]; 
ax.YLim = [0, mX_grand(2)]; 
ax.YTick = [];
ax.TickDir = 'out'; 

pnl.fontsize = 8;
pnl.margin = [5, 10, 5, 3]; 


if strcmp(noise_type, 'fractal')
    fname = sprintf('onlyNoise_noise-fractal_exp-%.1f_nrep-%d_mX', ...
                   noise_exponent, n_rep); 
elseif strcmp(noise_type, 'eeg')
    fname = sprintf('onlyNoise_noise-eeg_nrep-%d_mX', ...
                   n_rep); 
else
    error('noise type "%s" not implemented', noise_type);
end


if par.save_figs
   save(fullfile(par.data_path, fname), 'mX', 'freq'); 
   save_fig(f, fname)
end

   





























