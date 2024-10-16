clear all

% add the necessary library folders to MATLAB path (make 
% sure you change these to your local paths)
addpath(genpath('~/projects_git/rnb_tools'));
addpath(genpath('~/projects_git/acf_tools'));

% small utilities here 
addpath(genpath('lib')); 

%% 

n_rep = 100; 

n_events_pattern_all = 12 * [1, 2, 5, 10, 20]; 


f_fft = figure('Color', 'white'); 
pnl_fft = panel(f_fft); 
pnl_fft.pack('v', length(n_events_pattern_all))


f_acf = figure('Color', 'white'); 
pnl_acf = panel(f_acf); 
pnl_acf.pack('v', length(n_events_pattern_all))

xmax_fft = -inf; 
xmax_acf = -inf; 

colors = brewermap(10, 'Purples'); 

for i_n_events=1:length(n_events_pattern_all)
    
    
    var_fft = nan(1, n_rep); 
    var_acf = nan(1, n_rep); 

    for i_rep=1:n_rep


        % samping rate 
        fs = 200; 

        n_events_pattern = n_events_pattern_all(i_n_events); 

        n_events_total = 12 * 120; 


        % symbolic pattern of events on a fast isochronous time grid 
        pat = randsample([0, 1], n_events_pattern, true); 
        n_cycles = n_events_total / n_events_pattern; 

        % inter-onser interval between two grid points 
        grid_ioi = 0.2; 

        duty_cycle = 0.1; 

        ir = get_square_kernel(fs, ...
            'duration', duty_cycle, ...
            'rampon', 0, ...
            'rampoff', 0 ...
            ); 

        % make whole signal 
        [x, t] = get_s(...
                        pat, ...
                        grid_ioi, ...
                        fs, ...
                        'n_cycles', n_cycles, ...
                        'ir', ir, ...
                        'emph_magn', 0, ...
                        'emph_period', 4, ...
                        'jitter', 0 ...
                        );

        % frequencies of interst
        max_freq = 20; 
        f0_to_excl = 5; 

        [freq_meter_rel, freq_meter_unrel, frex] = get_meter_freq(max_freq, ...
                                                        'f0_to_excl', f0_to_excl);

        % lags of interest 
        signal_dur = length(x) / fs; 
        max_lag = signal_dur / 2; 

        lag_base_incl_meter_rel = [0.8]; 
        lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

        lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
        lag_base_excl_meter_unrel = [0.8]; 

        [lags_meter_rel, lags_meter_unrel] = get_meter_lags(...
                    max_lag, ...
                    lag_base_incl_meter_rel, lag_base_excl_meter_rel, ...
                    lag_base_incl_meter_unrel, lag_base_excl_meter_unrel ...
                    );

        % get acf withuout aperiodic subtraction    
        [acf, lags, ~, mX, freq] = get_acf(x, fs, 'normalize_x', true); 


        feat_fft = get_fft_features(mX, freq, ...
                                    freq_meter_rel, freq_meter_unrel); 

        feat_acf = get_acf_features(acf, lags, ...
                                    lags_meter_rel, lags_meter_unrel);    

        var_fft(i_rep) = var(feat_fft.vals); 
        var_acf(i_rep) = var(feat_acf.vals); 

    end
    
    

    %%

    fprintf('\n\n\npattern made from %d events\n',n_events_pattern); 


    xmax_fft = max([xmax_fft, var_fft]); 
    xmax_acf = max([xmax_acf, var_acf]); 
    
    ax = pnl_fft(i_n_events).select(); 
    hold(ax, 'on'); 
    h = histogram(ax, var_fft, 'DisplayStyle', 'bar', ...
                  'linew', 2, 'EdgeColor', colors(3+i_n_events, :), 'FaceColor', colors(3+i_n_events, :));
%     ax.XLim = [0, inf]; 
%     title(sprintf('var of fft at (N=%d) freqs', length(feat_fft.vals))); 

    ax = pnl_acf(i_n_events).select(); 
    hold(ax, 'on');  
    h = histogram(ax, var_acf, 'DisplayStyle', 'bar', ...
                  'linew', 2, 'EdgeColor', colors(3+i_n_events, :), 'FaceColor', colors(3+i_n_events, :));
%     ax.XLim = [0, inf]; 
%     title(sprintf('var of ACF at (N=%d) lags', length(feat_acf.vals))); 


end

                        
                        
                        
                        
for i_n_events=1:length(n_events_pattern_all)

    ax = pnl_fft(i_n_events).select(); 
    ax.XLim = [0, xmax_fft]; 
    
    ax = pnl_acf(i_n_events).select(); 
    ax.XLim = [0, xmax_acf]; 
end







