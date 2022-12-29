
clear 
fprintf('\n');
addpath(genpath('lib'))
addpath('/datadisk/Dropbox/projects/autocorrelation/acf_tools/src'); 

%% simulate

fs = 128; 
dur = 500; 
N = round(dur * fs); 

exponent = -1.5; 

snr = 0.2; 

fit_knee = false; 

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

%% 

ir = get_square_kernel(fs, ...
    'duration', 0.100, ...
    'rampon', 0, ...
    'rampoff', 0 ...
    ); 


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


%% repeat simulation - test bias/variance 

% no noise (clean signal) 
[acf, lags, ap, mX, freq] = get_acf(x_clean, fs); 

feat_acf_orig = get_acf_features(acf, lags, ...
                             lags_meter_rel, lags_meter_unrel, ...
                             'lags_meter_unrel_left', lags_meter_unrel_left, ...
                             'lags_meter_unrel_right', lags_meter_unrel_right); 

feat_fft_orig = get_fft_features(mX, freq, ...
                                 freq_meter_rel, freq_meter_unrel);
                             
    
% generate noisy signal (simulataneously for all repetitions)
noise = get_colored_noise2([n_rep, length(x_clean)], fs, exponent); 

x = x_clean + noise./snr; 

    
% withuout aperiodic subtraction    
[acf, lags, ap, mX, freq] = get_acf(x, fs); 
    
feat_acf = get_acf_features(acf, lags, ...
                            lags_meter_rel, lags_meter_unrel, ...
                            'lags_meter_unrel_left', lags_meter_unrel_left, ...
                            'lags_meter_unrel_right', lags_meter_unrel_right);

feat_fft = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
            
            
% with aperiodic subtraction    
acf_subtracted = get_acf(x, fs, 'rm_ap', true,  'f0_to_ignore', 1/2.4); 

feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                 lags_meter_rel, lags_meter_unrel, ...
                                 'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                 'lags_meter_unrel_right', lags_meter_unrel_right); 
                                   
mX_subtracted = subtract_noise_bins(mX, noise_bins(1),  noise_bins(2)); 

feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                       freq_meter_rel, freq_meter_unrel);
                                       

% plot 
% ----

col_raw = [152, 63, 212]/255;
col_subtracted = [17, 120, 48]/255;


f = figure('color', 'white', 'position', [945 545 430 199], 'name', 'faetures'); 
pnl = panel(f); 

pnl.pack('h', 4); 

ax = pnl(1).select(); 

feat = []; 
feat(1).name = 'raw'; 
feat(1).data = feat_acf.ratio_meter_rel; 
feat(1).color = col_raw; 

feat(2).name = 'subtracted'; 
feat(2).data = feat_acf_subtracted.ratio_meter_rel; 
feat(2).color = col_subtracted; 

plot_multiple_cond('ax', ax, ...
                  'plot_legend', true, ...
                  'feat', feat, ...
                  'feat_orig', feat_acf_orig.ratio_meter_rel); 

pnl(1).ylabel('ACF ratio'); 


ax = pnl(2).select(); 

feat = []; 
feat(1).name = 'raw'; 
feat(1).data = feat_acf.ratio_meter_rel_left; 
feat(1).color = col_raw; 

feat(2).name = 'subtracted'; 
feat(2).data = feat_acf_subtracted.ratio_meter_rel_left; 
feat(2).color = col_subtracted; 

plot_multiple_cond('ax', ax, ...
                  'plot_legend', false, ...
                  'feat', feat, ...
                  'feat_orig', feat_acf_orig.ratio_meter_rel_left); 

pnl(2).ylabel('ACF ratio L'); 


ax = pnl(3).select(); 

feat = []; 
feat(1).name = 'raw'; 
feat(1).data = feat_acf.ratio_meter_rel_right; 
feat(1).color = col_raw; 

feat(2).name = 'subtracted'; 
feat(2).data = feat_acf_subtracted.ratio_meter_rel_right; 
feat(2).color = col_subtracted; 

plot_multiple_cond('ax', ax, ...
                  'plot_legend', false, ...
                  'feat', feat, ...
                  'feat_orig', feat_acf_orig.ratio_meter_rel_right); 

pnl(3).ylabel('ACF ratio R'); 


ax = pnl(4).select(); 

feat = []; 
feat(1).name = 'raw'; 
feat(1).data = feat_fft.z_meter_rel; 
feat(1).color = col_raw; 

feat(2).name = 'subtracted'; 
feat(2).data = feat_fft_subtracted.z_meter_rel; 
feat(2).color = col_subtracted; 

plot_multiple_cond('ax', ax, ...
                  'plot_legend', false, ...
                  'feat', feat, ...
                  'feat_orig', feat_fft_orig.z_meter_rel); 

pnl(4).ylabel('FFT z meter'); 

pnl(4).marginleft = 25; 
pnl.margin = [15, 5, 5, 15]; 






