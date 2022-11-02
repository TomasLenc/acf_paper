function f = plot_example(x, t, acf, lags, ap, mX, freq, ...
                      lags_meter_rel, lags_meter_unrel, ...
                      freq_meter_rel, freq_meter_unrel, varargin) 

parser = inputParser; 

addParameter(parser, 'lags_meter_unrel_left', []); 
addParameter(parser, 'lags_meter_unrel_right', []); 
addParameter(parser, 'mX_subtr', []); 
addParameter(parser, 'acf_subtr', []); 
addParameter(parser, 'max_lag', 1.2); 
addParameter(parser, 'max_freq', 5.2); 
addParameter(parser, 'time_col', [0, 0, 0]); 
addParameter(parser, 'prec', 1000); 

parse(parser, varargin{:}); 

mX_subtr = parser.Results.mX_subtr; 
acf_subtr = parser.Results.acf_subtr; 

lags_meter_unrel_left = parser.Results.lags_meter_unrel_left; 
lags_meter_unrel_right = parser.Results.lags_meter_unrel_right; 

prec = parser.Results.prec; 
max_lag = parser.Results.max_lag; 
max_freq = parser.Results.max_freq; 

time_col = parser.Results.time_col; 


%%
                  
feat_acf = get_acf_features(acf, lags,...
                            lags_meter_rel, lags_meter_unrel, ...
                            'lags_meter_unrel_left', lags_meter_unrel_left, ...
                            'lags_meter_unrel_right', lags_meter_unrel_right); 
                        
feat_fft = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel);

if ~isempty(acf_subtr)
    feat_acf_subtr = get_acf_features(acf_subtr, lags,...
                                lags_meter_rel, lags_meter_unrel, ...
                                'lags_meter_unrel_left', lags_meter_unrel_left, ...
                                'lags_meter_unrel_right', lags_meter_unrel_right); 
end

if ~isempty(mX_subtr)
    feat_fft_subtracted = get_fft_features(mX_subtr, freq, ...
                                    freq_meter_rel, freq_meter_unrel);
end


%%

% open figure
f = figure('color','white', 'position', [673 485 1062 240]); 

pnl = panel(f); 

pnl.pack('h', [50, 25, 25]); 
pnl(1).pack({[0, 0, 1, 1]}); 
pnl(2).pack({[0, 0, 1, 1]}); 
pnl(3).pack({[0, 0, 1, 1]}); 

inset_coord = [0.5, 1.5, 0.5, 0.65]; 
pnl(1).pack({inset_coord});
pnl(2).pack({inset_coord});
pnl(3).pack({inset_coord});


% plot time-domain 
ax = pnl(1, 1).select(); 
plot_time(ax, x, t, 'col', time_col)
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
         'max_freq', max_freq); 
ax.YTick = []; 

if ~isempty(mX_subtr)
    ax = pnl(2, 2).select(); 
    features = []; 
    features.z = feat_fft_subtracted.z_meter_rel; 
    plot_fft(ax, mX_subtr, freq, ...
             'freq_meter_rel', freq_meter_rel, ...
             'freq_meter_unrel', freq_meter_unrel, ...
             'features', features, ...
             'max_freq', max_freq); 
    ax.XAxis.Visible = 'off';  
    ax.YAxis.Visible = 'off';  
end


% plot ACF
ax = pnl(3, 1).select(); 

features = []; 
features.ratio = feat_acf.ratio_meter_rel; 
features.contr = feat_acf.contrast_meter_rel; 

idx = dsearchn(lags', max_lag); 
acf_to_plot = (acf - min(acf(1:idx))) ./ (max(acf(1:idx)) - min(acf(1:idx))); 

plot_acf(ax, ...
         acf_to_plot, ...
         lags, ...
         'features', features, ...
         'lags_meter_rel', lags_meter_rel, ...
         'lags_meter_unrel', lags_meter_unrel, ...
         'max_lag', max_lag, ...
         'prec', prec); 
ax.YTick = []; 

if ~isempty(acf_subtr)
    ax = pnl(3, 2).select(); 
    features = []; 
    features.ratio = feat_acf_subtr.ratio_meter_rel; 
    features.contr = feat_acf_subtr.contrast_meter_rel; 
    acf_to_plot = (acf_subtr - min(acf_subtr(1:idx))) ./ ...
                  (max(acf_subtr(1:idx)) - min(acf_subtr(1:idx))); 
    plot_acf(ax, ...
             acf_to_plot, ...
             lags, ...
             'features', features, ...
             'lags_meter_rel', lags_meter_rel, ...
             'lags_meter_unrel', lags_meter_unrel, ...
             'max_lag', max_lag, ...
             'prec', prec); 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
end


% margins / labels
pnl(1).xlabel('time (s)')
pnl(2).xlabel('frequency (Hz)')
pnl(3).xlabel('lag (s)')

pnl.de.margin = [15, 10, 10, 15]; 
pnl(2).marginleft = 9; 
pnl(3).marginleft = 20; 
pnl.margin = [15, 12, 25, 31]; 

pnl.fontsize = 12;
