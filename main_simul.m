% creates rhythmic signals to test autocorrelation method
clear all

addpath('lib')

plot_fft = false; 


%% set signal parameters
% each parameter is set for both signals that will be compared

track_name = {'stim1', 'stim2'}; 
   

cols = {[152, 63, 212]/255, 
        [17, 120, 48]/255};

% set rhythmic patterns  
pattern = {[1 0 1 1 1 1 0 1 1 1 0 0]
           [1 0 1 1 1 1 0 1 1 1 0 0]};

% sampling rate
fs = 1000;

% time interval between successive events (either sound or silence)
grid_ioi = 0.2; 

% how many times the rhythmic pattern repeats in each trial
n_cycles = 16; % 40.8 sec is 17 cycles (for 2.4s cycle)

emph_magn = [0, 0]; 

emph_period = [4, 4]; 

jitter = [0, 0] * grid_ioi; 

snr = [Inf, Inf]; % Inf, 0.7

noise_exp = [-2, -2]; 



% number of simulated repetitions 
n_rep = 1; 

% whether to save continuous waveforms (s, acf, ...) for plotting 
% set to false if simulating too many repetitions
save_continuous = true; 

% whether to subtract the estimated aperiodic component from the acf before
% calculating indices of meter periodicities
subtract_ap = false; 

% whether to fit the apeirodic component with a knee
fit_knee = false; 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
lags_meter_rel = [0.8]; 
lags_meter_unrel = [0.6, 1.0]; 

% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 

max_lag = 1.2; 

% =========================================================================
%%%%% set impulse response parameters
ir_par = {}; 

% %%% square
% tmp = []; 
% % type of impulse response (provide parameters accordingly)
% tmp.type = 'square'; 
% % duration of sound event
% tmp.eventDur = 0.200; % 4 ms for click, 100ms for tone
% % duration of linear onset ramp for the sound event
% tmp.rampon = 0.100;
% % duration of linear offset ramp for the sound event
% tmp.rampoff = 0.100; 

% %%% click
% tmp = []; 
% % type of impulse response (provide parameters accordingly)
% tmp.type = 'click'; 
% % duration of sound event
% tmp.eventDur     = 0.004; % 4 ms for click, 100ms for tone

% %%% erp
% tmp = struct('type','erp', ...
%                 'A',[0.4,0.75], ...
%                 'T0',[0,0], ...
%                 'Tau',[0.2, 0.050], ...
%                 'F0',[1,7], ...
%                 'Dur',[0.5]);
% 

% --------- IR1 ---------

%%% square
tmp = []; 
% type of impulse response (provide parameters accordingly)
tmp.type = 'square'; 
% duration of sound event
tmp.eventDur = 0.050; % 4 ms for click, 100ms for tone
% duration of linear onset ramp for the sound event
tmp.rampon = 0.0;
% duration of linear offset ramp for the sound event
tmp.rampoff = 0.0; 

ir_par{1} = tmp; 

% --------- IR2 ---------

%%% square
tmp = []; 
% type of impulse response (provide parameters accordingly)
tmp.type = 'square'; 
% duration of sound event
tmp.eventDur = 0.180; % 4 ms for click, 100ms for tone
% duration of linear onset ramp for the sound event
tmp.rampon = 0.0;
% duration of linear offset ramp for the sound event
tmp.rampoff = 0.0; 

ir_par{2} = tmp; 



%% synthesis + analysis 

res = []; 
% go over tracks
for i_cond=1:length(track_name)

    % create impulse response
    ir = get_ir(ir_par{i_cond}, fs); 
    
    % make whole signal 
    [s_clean, t] = get_s(...
                        pattern{i_cond}, ...
                        n_cycles, ...
                        grid_ioi, ...
                        jitter(i_cond), ...
                        ir, ...
                        fs, ...
                        'emph_magn', emph_magn(i_cond), ...
                        'emph_period', emph_period(i_cond) ...
                        );
                    
    % save params
    res(i_cond).track_name  = track_name{i_cond}; 
    res(i_cond).pattern     = pattern{i_cond}; 
    res(i_cond).grid_ioi    = grid_ioi; 
    res(i_cond).jitter      = jitter(i_cond); 
    res(i_cond).ir_par      = ir_par{i_cond}; 
    res(i_cond).n_cycles    = n_cycles; 
    res(i_cond).trial_dur   = length(s_clean) / fs; 
    res(i_cond).fs          = fs;
    res(i_cond).t           = t;
    
    cycle_dur = length(res(i_cond).pattern) * res(i_cond).grid_ioi; 

    if save_continuous
        res(i_cond).s           = nan(n_rep, length(s_clean));
        res(i_cond).acf_raw_    = nan(n_rep, length(s_clean) / 2);
        res(i_cond).acf_        = nan(n_rep, length(s_clean) / 2); 
        res(i_cond).lags_time_  = nan(n_rep, length(s_clean) / 2); 
        res(i_cond).ap_         = nan(n_rep, length(s_clean) / 2); 
    end
    res(i_cond).ap_par_                     = cell(1, n_rep); 
    res(i_cond).z_meter_rel_                = nan(1, n_rep); 
    res(i_cond).ratio_meter_rel_            = nan(1, n_rep); 
    res(i_cond).contrast_meter_rel_         = nan(1, n_rep); 
    res(i_cond).contrast_meter_rel_left_    = nan(1, n_rep); 
    res(i_cond).contrast_meter_rel_right_   = nan(1, n_rep); 

    for i_rep=1:n_rep

        fprintf('%s: %d/%d\n', track_name{i_cond}, i_rep, n_rep); 
        
        s = add_scaled_noise(s_clean, snr(i_cond), fs, noise_exp(i_cond)); 
        
        [acf_raw, lags_time] = get_acf(s, res(i_cond).fs);
        
        [acf, ~, ap, ~, ~, ap_par] = get_acf(s, res(i_cond).fs, ...
                               'rm_ap', true, ...
                               'f0_to_ignore', 1/2.4);

        % get indices for lags of interest
        lags_meter_rel_idx = dsearchn(lags_time', lags_meter_rel'); 
        lags_meter_unrel_idx = dsearchn(lags_time', lags_meter_unrel'); 
        lags_meter_unrel_left_idx = dsearchn(lags_time', lags_meter_unrel_left'); 
        lags_meter_unrel_right_idx = dsearchn(lags_time', lags_meter_unrel_right'); 

        % calculate mean acf value for lags of interest
        acf_mean_meter_rel = mean(acf(lags_meter_rel_idx)); 
        acf_mean_meter_unrel = mean(acf(lags_meter_unrel_idx)); 
        acf_mean_meter_unrel_left = mean(acf(lags_meter_unrel_left_idx)); 
        acf_mean_meter_unrel_right = mean(acf(lags_meter_unrel_right_idx)); 

        % z-score
        z = zscore([acf(lags_meter_rel_idx), acf(lags_meter_unrel_idx)]); 
        z_meter_rel = mean(z(1:length(acf_mean_meter_rel))); 

        % acf ratio 
        ratio_meter_rel = acf_mean_meter_rel / acf_mean_meter_unrel;
        ratio_meter_rel_left = acf_mean_meter_rel / acf_mean_meter_unrel_left;
        ratio_meter_rel_right = acf_mean_meter_rel / acf_mean_meter_unrel_right;

        % acf contrast 
        contrast_meter_rel = (acf_mean_meter_rel-acf_mean_meter_unrel) / ...
                             (acf_mean_meter_rel+acf_mean_meter_unrel); 

        % save results into structure
        if save_continuous
            res(i_cond).s(i_rep, :)         = s;
            res(i_cond).acf_raw_(i_rep, :)  = acf_raw; 
            res(i_cond).acf_(i_rep, :)      = acf; 
            res(i_cond).ap_(i_rep, :)       = ap; 
        end
        
        res(i_cond).ap_par_(i_rep)                  = ap_par; 
        res(i_cond).z_meter_rel_(i_rep)             = z_meter_rel; 
        res(i_cond).ratio_meter_rel_(i_rep)         = ratio_meter_rel; 
        res(i_cond).ratio_meter_rel_left_(i_rep)    = ratio_meter_rel_left; 
        res(i_cond).ratio_meter_rel_right_(i_rep)   = ratio_meter_rel_right; 
        res(i_cond).contrast_meter_rel_(i_rep)      = contrast_meter_rel; 
        
    end
    
    res(i_cond).lags_time_ = lags_time; 
    
end

%% 

if plot_fft

    f = figure('color', 'white', 'name', 'FFT'); 
    pnl = panel(f); 
    pnl.pack('v', 2); 

    for i_cond=1:length(track_name)

        powX = pwelch(res(1).s(1, :)); 
        powX = powX .* 0; 
        for i_rep=1:n_rep
            [powX_tmp, freq_norm] = pwelch(res(i_cond).s(i_rep, :)); 
            powX = [powX + powX_tmp]; 
        end
        powX = powX ./ n_rep; 

        mX = sqrt(powX); 
        freq = freq_norm / (2*pi) * res(i_cond).fs; 

        pnl(1).select(); 
        plot(freq, mX, 'linew', 2, 'color', cols{i_cond})
        hold on
        title('magnitude')
        xlim([0.3, 30])

        pnl(2).select(); 
        plot(freq, powX, 'linew', 2, 'color', cols{i_cond})
        hold on
        title('power')
        xlim([0.3, 30])

    end
    
end




%% plot ratios

if length(res) == 2 && n_rep > 1
    
    f = figure('color', 'white', 'position', [945 545 430 199], 'name', 'faetures'); 
    pnl = panel(f); 

    pnl.pack('h', 3); 

    ax = pnl(1).select(); 
    plot_points_cond(ax, ...
                     res(1).ratio_meter_rel_, ...
                     res(2).ratio_meter_rel_, ...
                     cols{1}, ...
                     cols{2})
    ylabel(ax, 'ratio_meter_rel', 'Interpreter', 'none'); 


    ax = pnl(2).select(); 

    plot_points_cond(ax, ...
                     res(1).ratio_meter_rel_left_, ...
                     res(2).ratio_meter_rel_left_, ...
                     cols{1}, ...
                     cols{2})
    ylabel(ax, 'ratio_meter_rel_left', 'Interpreter', 'none'); 


    ax = pnl(3).select(); 
    plot_points_cond(ax, ...
                     res(1).ratio_meter_rel_right_, ...
                     res(2).ratio_meter_rel_right_, ...
                     cols{1}, ...
                     cols{2})
    ylabel(ax, 'ratio_meter_rel_right', 'Interpreter', 'none'); 


    pnl(2).marginleft = 30; 
    pnl.margin = [20, 10, 5, 5]; 


end


%% plot continuous waveforms 

% set some parameters 
prec = 1e6; 

col_acf = repmat(0, 1, 3);  
col_ap = [224, 117, 29]/255;  
col_meter_rel = [214, 52, 24]/255; 
col_meter_unrel = [45, 114, 224]/255; 

f = figure('color','white', 'position', [673, 525, 1062, 200 * length(res)], ...
           'name', 'acf'); 
pnl = panel(f); 

pnl.pack('h', 3)
pnl(1).pack('v', length(res))
pnl(2).pack('v', length(res))
pnl(3).pack('v', length(res))

for i_cond=1:length(res)

    % plot time-domain 
    ax = pnl(1, i_cond).select(); 
    plot_time(ax, res(i_cond).s(1, :), res(i_cond).t, 'col', cols{i_cond}); 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 

    % plot autocorrelation function
    ax = pnl(2, i_cond).select(); 
    plot_acf(ax, ...
             res(i_cond).acf_raw_(1, :), ...
             res(i_cond).lags_time_, ...
             'ap', res(i_cond).ap_(1, :), ...
             'lags_meter_rel', lags_meter_rel, ...
             'lags_meter_unrel', lags_meter_unrel, ...
             'max_lag', max_lag, ...
             'col_acf', col_acf, ...
             'col_ap', col_ap, ...
             'col_meter_rel', col_meter_rel, ...
             'col_meter_unrel', col_meter_unrel, ...
             'prec', prec); 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 

    % plot corrected autocorrelation function
    ax = pnl(3, i_cond).select(); 
    features = []; 
    features.z = res(i_cond).z_meter_rel_(1); 
    features.ratio = res(i_cond).ratio_meter_rel_(1); 
    features.contr = res(i_cond).contrast_meter_rel_(1); 
    plot_acf(ax, ...
             res(i_cond).acf_(1, :), ...
             res(i_cond).lags_time_, ...
             'features', features, ...
             'lags_meter_rel', lags_meter_rel, ...
             'lags_meter_unrel', lags_meter_unrel, ...
             'max_lag', max_lag, ...
             'col_acf', col_acf, ...
             'col_ap', col_ap, ...
             'col_meter_rel', col_meter_rel, ...
             'col_meter_unrel', col_meter_unrel, ...
             'prec', prec); 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
      
end

pnl(1).xlabel('time (s)')
pnl(1).ylabel('amplitude')

pnl(2).xlabel('lag (s)')
pnl(2).ylabel('autocorrelation')

pnl(3).xlabel('lag (s)')
pnl(3).ylabel('autocorrelation')

pnl.fontsize = 12; 


%% save figure

% savePath = 'output/acf_syncopated_squareVsTriangle_noies'; 
% 
% saveas(f, sprintf([savePath,'.fig'])) 
% print(f, '-depsc', '-painters', [savePath,'.eps']); 

