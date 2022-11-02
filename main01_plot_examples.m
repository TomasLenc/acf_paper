
clear 
fprintf('\n');
addpath(genpath('lib'))
addpath('/datadisk/Dropbox/projects/autocorrelation/acf_tools/src'); 

%% simulate

fs = 200; 

exponent = -1.4; 

snr = 1; 

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

ir_par = []; 
ir_par.type = 'square'; 
ir_par.eventDur = 0.030;
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

x_clean = x_clean - min(x_clean); 

noise = get_colored_noise(length(x_clean), fs, exponent); 
noise = (noise - mean(noise)) ./ std(noise); 
x = x_clean + noise/snr; 

x = x - min(x); 

%% 

for i_cond=1:3
    
    switch i_cond
        
        case 1
            % no noise (clean signal) 
            [acf, lags, ap, mX, freq] = get_acf(x_clean, fs); 
            x_to_plot = x_clean; 

        case 2
            % withuout aperiodic subtraction
            [acf, lags, ap, mX, freq] = get_acf(x, fs); 
            x_to_plot = x; 
            
        case 3
            % with aperiodic subtraction
            [acf, lags, ap, mX, freq] = get_acf(x, fs, ...
                                               'rm_ap', true, ...
                                               'f0_to_ignore', 1/2.4 ...
                                               ); 
            x_to_plot = x; 
            
    end
            
%     % normalize acf
%     max_idx = dsearchn(lags', 1.2); 
%     lags = lags(1:max_idx); 
%     acf = acf(1:max_idx); 
%     acf_min = min(acf); 
%     acf_max = max(acf); 
%     acf_range = acf_max - acf_min; 
%     acf = (acf - acf_min) ./ acf_range; 

    feat_acf = get_acf_features(acf, lags, ...
                                lags_meter_rel, lags_meter_unrel, ...
                                'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                'lags_meter_unrel_right', lags_meter_unrel_right); 


    mX_subtracted = subtract_noise_bins(mX, noise_bins(1),  noise_bins(2)); 

    feat_fft = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel);

    feat_fft_subtracted = get_fft_features(mX_subtracted, freq, freq_meter_rel, freq_meter_unrel);


    % open figure
    f = figure('color','white', 'position', [673 485 1062 240]); 

    pnl = panel(f); 

    pnl.pack('h', [50, 25, 25]); 
    pnl(2).pack({[0, 0, 1, 1]}); 

    inset_coord = [0.7, 1.0, 0.4, 0.65]; 
    pnl(2).pack({inset_coord});


    % plot time-domain 
    ax = pnl(1).select(); 
    plot_time(ax, x_to_plot, t)
    ax.YAxis.Visible = 'off'; 


    % plot FFT
    ax = pnl(2, 1).select(); 
    features = []; 
    features.z = feat_fft.z_meter_rel; 
    plot_fft(ax, mX, freq, ...
             'ap', ap, ...
             'freq_meter_rel', freq_meter_rel, ...
             'freq_meter_unrel', freq_meter_unrel, ...
             'features', features, ...
             'max_freq', max_freq_plot); 
    ax.YTick = []; 

    ax = pnl(2, 2).select(); 
    features = []; 
    features.z = feat_fft_subtracted.z_meter_rel; 
    plot_fft(ax, mX_subtracted, freq, ...
             'freq_meter_rel', freq_meter_rel, ...
             'freq_meter_unrel', freq_meter_unrel, ...
             'features', features, ...
             'max_freq', max_freq_plot); 
    ax.XAxis.Visible = 'off';  
    ax.YAxis.Visible = 'off';  


    % plot ACF
    ax = pnl(3).select(); 
    features = []; 
    features.ratio = feat_acf.ratio_meter_rel; 
    features.ratio_L = feat_acf.ratio_meter_rel_left; 
    features.ratio_R = feat_acf.ratio_meter_rel_right; 

    idx = dsearchn(lags', 1.2); 
    acf_to_plot = (acf - min(acf(1:idx))) ./ (max(acf(1:idx)) - min(acf(1:idx))); 
    
    plot_acf(ax, ...
             acf_to_plot, ...
             lags, ...
             'features', features, ...
             'lags_meter_rel', lags_meter_rel, ...
             'lags_meter_unrel', lags_meter_unrel, ...
             'max_lag', 1.2, ...
             'prec', 1000); 
    ax.YTick = []; 

    pnl(1).xlabel('time (s)')
    pnl(2).xlabel('frequency (Hz)')
    pnl(3).xlabel('lag (s)')

    pnl.de.margin = [15, 10, 10, 15]; 
    pnl(2).marginleft = 9; 
    pnl(3).marginleft = 20; 
    pnl.margin = [15, 15, 16, 25]; 

    pnl.fontsize = 12;




end



