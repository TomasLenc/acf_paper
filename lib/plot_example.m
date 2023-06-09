function f = plot_example(x, t, acf, lags, ap, mX, freq, ...
                      lags_meter_rel, lags_meter_unrel, ...
                      freq_meter_rel, freq_meter_unrel, varargin) 

parser = inputParser; 

addParameter(parser, 'pnl', []); 
addParameter(parser, 'mX_subtr', []); 
addParameter(parser, 'acf_subtr', []); 
addParameter(parser, 'normalize_acf_for_plotting', true); 
addParameter(parser, 'max_t', 4.8); 
addParameter(parser, 'min_lag', min(lags)); 
addParameter(parser, 'max_lag', 1.21); 
addParameter(parser, 'max_freq', 5.2); 
addParameter(parser, 'time_col', [0, 0, 0]); 
addParameter(parser, 'subplot_proportions', [50, 25, 25]); 
addParameter(parser, 'plot_time_xaxis', true); 
addParameter(parser, 'plot_xticks', true); 
addParameter(parser, 'plot_xlabels', true); 
addParameter(parser, 'plot_features', true); 
addParameter(parser, 'fontsize', 12); 
addParameter(parser, 'prec', 1000); 


parse(parser, varargin{:}); 

pnl = parser.Results.pnl; 

mX_subtr = parser.Results.mX_subtr; 
acf_subtr = parser.Results.acf_subtr; 
prec = parser.Results.prec; 
subplot_proportions = parser.Results.subplot_proportions; 
plot_time_xaxis = parser.Results.plot_time_xaxis; 
plot_xticks = parser.Results.plot_xticks; 
plot_xlabels = parser.Results.plot_xlabels; 
plot_features = parser.Results.plot_features; 
max_t = parser.Results.max_t; 
min_lag = parser.Results.min_lag; 
max_lag = parser.Results.max_lag; 
max_freq = parser.Results.max_freq; 
fontsize = parser.Results.fontsize; 
normalize_acf_for_plotting = parser.Results.normalize_acf_for_plotting; 

time_col = parser.Results.time_col; 


%%
                  
feat_acf = get_acf_features(acf, lags,...
                            lags_meter_rel, lags_meter_unrel); 
                        
feat_fft = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel);

if ~isempty(acf_subtr)
    feat_acf_subtr = get_acf_features(acf_subtr, lags,...
                                lags_meter_rel, lags_meter_unrel); 
end

if ~isempty(mX_subtr)
    feat_fft_subtracted = get_fft_features(mX_subtr, freq, ...
                                    freq_meter_rel, freq_meter_unrel);
end


%%

% open figure
if isempty(pnl)
    f = figure('color','white', 'position', [95 67 1062 240]); 
    pnl = panel(f); 
else
    f = []; 
    pnl.pack('h', subplot_proportions); 
    pnl(1).pack({[0, 0, 1, 1]}); 
    pnl(2).pack({[0, 0, 1, 1]}); 
    pnl(3).pack({[0, 0, 1, 1]}); 

    inset_coord = [0.3, 1.22, 0.7, 0.71]; 
    pnl(1).pack({inset_coord});
    pnl(2).pack({inset_coord});
    pnl(3).pack({inset_coord});
end


% plot time-domain 
ax = pnl(1, 1).select(); 
plot_erp(x, 't', t, 'col', time_col, 'ax', ax);
ax.XLim = [0, min(t(end), max_t)]; 
ax.YAxis.Visible = 'off'; 
if ~plot_time_xaxis
    ax.XAxis.Visible = 'off'; 
end
if ~plot_xticks
    ax.XTick = []; 
end


% plot raw FFT
ax = pnl(2, 1).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'frex_meter_rel', freq_meter_rel, ...
         'frex_meter_unrel', freq_meter_unrel, ...
         'maxfreqlim', max_freq); 
ax.YTick = []; 
if ~plot_xticks
    ax.XTick = []; 
end


% plot 1/f component
if ~isempty(ap)
    hold(ax, 'on');
    plot(ax, freq, ap, '--', 'color', 'k', 'linew', 2);    
end


% plot calculated feature values 
if plot_features
    features = []; 
    features.z = feat_fft.z_meter_rel; 

    tit = ''; 
    keys = fieldnames(features); 
    for i_key=1:length(keys)
        tit = [tit, sprintf('%s=%.2f   ', keys{i_key}, features.(keys{i_key}))];    
    end
    h_tit = title(ax, tit, 'Interpreter', 'none');
    h_tit.HorizontalAlignment = 'left'; 
    h_tit.Units = 'normalized'; 
    h_tit.Position(1) = 0; 
    h_tit.Position(2) = 1.1; 
end

% plot SNR-subtracted FFT
if ~isempty(mX_subtr)
    ax = pnl(2, 2).select(); 
    features = []; 
    features.z = feat_fft_subtracted.z_meter_rel; 
    plot_fft(freq, mX_subtr, ...
             'ax', ax, ...
             'frex_meter_rel', freq_meter_rel, ...
             'frex_meter_unrel', freq_meter_unrel, ...
             'maxfreqlim', max_freq); 
    ax.XAxis.Visible = 'off';  
    ax.YAxis.Visible = 'off';  
end


% plot ACF
ax = pnl(3, 1).select(); 

if max_lag > 10
    linew_acf = 0.8;
else
    linew_acf = 2;
end

features = struct;
if plot_features
    features.ratio = feat_acf.ratio_meter_rel; 
    features.z = feat_acf.z_meter_rel; 
end

if normalize_acf_for_plotting
    idx = dsearchn(lags', max_lag); 
    acf_to_plot = (acf - min(acf(1:idx))) ./ (max(acf(1:idx)) - min(acf(1:idx))); 
else
    acf_to_plot = acf; 
end

plot_acf(ax, ...
         acf_to_plot, ...
         lags, ...
         'features', features, ...
         'lags_meter_rel', lags_meter_rel, ...
         'lags_meter_unrel', lags_meter_unrel, ...
         'min_lag', min_lag, ...
         'max_lag', max_lag, ...
         'linew', linew_acf, ...
         'opacity_lagz', 0.5, ...
         'prec', prec); 
ax.YTick = []; 
if ~plot_xticks
    ax.XTick = []; 
end

if ~isempty(acf_subtr)
    ax = pnl(3, 2).select(); 
    features = struct; 
    if plot_features
        features.ratio = feat_acf_subtr.ratio_meter_rel; 
        features.z = feat_acf_subtr.z_meter_rel; 
    end
    if normalize_acf_for_plotting
        idx = dsearchn(lags', max_lag); 
        acf_to_plot = (acf_subtr - min(acf_subtr(1:idx))) ./ ...
                      (max(acf_subtr(1:idx)) - min(acf_subtr(1:idx))); 
    else
        acf_to_plot = acf_subtr; 
    end
    plot_acf(ax, ...
             acf_to_plot, ...
             lags, ...
             'features', features, ...
             'lags_meter_rel', lags_meter_rel, ...
             'lags_meter_unrel', lags_meter_unrel, ...
             'min_lag', min_lag, ...
             'max_lag', max_lag, ...
             'linew', linew_acf, ...
             'opacity_lagz', 0.5, ...
             'prec', prec); 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
end


% margins / labels
if plot_xlabels
    pnl(1).xlabel('time (s)')
    pnl(2).xlabel('frequency (Hz)')
    pnl(3).xlabel('lag (s)')
end

pnl.de.margin = [15, 10, 10, 15]; 
pnl(2).marginleft = 9; 
pnl(3).marginleft = 9; 
pnl.margin = [15, 12, 10, 36]; 

pnl.fontsize = fontsize;
