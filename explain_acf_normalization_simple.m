addpath(genpath('/datadisk/projects_git_dl/rnb_tools/src'));
addpath(genpath('/datadisk/projects_git_dl/acf_tools/src'));
addpath(genpath('lib'))


%%


sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

% frequencies of interst
par.max_freq = 5; 
par.max_freq_plot = 5.1; 
par.f0_to_excl = 5; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq, ...
                                                'f0_to_excl', par.f0_to_excl);

% lags of interest 
par.max_lag = par.trial_dur / 2; 

par.lag_base_incl_meter_rel = [0.8]; 
par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
par.lag_base_excl_meter_unrel = [0.4]; 

[par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
            par.max_lag, ...
            par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
            par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
            );
        



%%
pat = [1 1 1 1 0 1 1 1 0 0 1 0];

grid_ioi = 0.2;

fs = 100;

n_cycles = 16;

ir = get_square_kernel(fs);

noise_exponent = -1.5; 

snr_target = 0.7;

response_f0 = 1 / (length(pat) * grid_ioi);

%%

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    pat, ...
                    grid_ioi, ...
                    fs, ...
                    'n_cycles', n_cycles, ...
                    'ir', ir ...
                    );
                
noise = get_colored_noise(length(x_clean), fs, noise_exponent); 
        
% scale the noise to the correct SNR 
x = add_signal_noise(x_clean, noise, snr_target);

%% estimate 1/f

% get whole-trial FFT
N = size(x, ndims(x)); 
hN = floor(N / 2) + 1; 
freq = [0 : hN-1]' / N * fs; 
mX = abs(fft(x, [], ndims(x))) / N * 2; 
mX = mX(1:hN); 

% get lags for ACF
lags = [0 : hN-1] / fs; 

% find frequencies which will be considered for the 1/f fit
min_freq_idx = dsearchn(freq, 0.1); 
max_freq_idx = dsearchn(freq, 9); 
freq_to_fit = freq(min_freq_idx : max_freq_idx); 

% ignore all harmonics of f0 up to nyquist frequency
nyq = fs/2; 
freq_to_ignore = [response_f0 : response_f0 : nyq]'; 
freq_to_ignore_idx = dsearchn(freq, freq_to_ignore); 

% for 1/f fitting, replace harmonics of f0 with mean of the bins around
bins = [2, 5];

mX_to_fit = mX; 
for i_f=1:length(freq_to_ignore)
    idx_1 = max(freq_to_ignore_idx(i_f) - bins(2), 1); 
    idx_2 = max(freq_to_ignore_idx(i_f) - bins(1), 1); 
    idx_3 = min(freq_to_ignore_idx(i_f) + bins(1), hN); 
    idx_4 = min(freq_to_ignore_idx(i_f) + bins(2), hN); 
    index = [idx_1:idx_2, idx_3:idx_4];
    
    mean_around_bin = mean(mX(index)); 
    mX_to_fit(freq_to_ignore_idx(i_f)) = mean_around_bin; 
end

% get log power spectra
log_pow = log10(mX_to_fit .^ 2); 

% restrict only to frequency range for 1/f fitting
log_pow_to_fit = log_pow(min_freq_idx : max_freq_idx); 

% fit aperiodic component
log_pow_to_fit = ensure_col(log_pow_to_fit);

[ap_par, ap_optim_exitflag] = fit_aperiodic(...
                    freq_to_fit, ...
                    log_pow_to_fit, ...
                    'robust', false, ...
                    'fit_knee', false,...
                    'verbose', true ...
                    ); 

% generate aperiodic with the estimated parameters across the whole
% frequency range 
ap = aperiodic(ap_par, freq)'; 

% convert from loq-power space to linear magnitude space 
ap_pow = 10.^ap; 
ap_linear = sqrt(ap_pow); 


%% remove 1/f and estimate ACF


ap_linear(isinf(ap_linear)) = 0; 

ap_for_norm = ap_linear; 

% mirror the aperiodic compoent so we also have it for negative frequencies 
if mod(N, 2) == 0
    % even number of frequency bins => we have bin at pi/2 => make sure you
    % don't copy this bin twice! 
    index = [2 : (size(ap_for_norm, ndims(ap_for_norm)) - 1)];      
else
    % odd number of frequency bins => no bin at pi/2 => just skip DC and mirror
    index = [2 : size(ap_for_norm, ndims(ap_for_norm))];      
end
ap_whole_spect = cat(...
                     ndims(x), ...
                     ap_for_norm, ...
                     flip(ap_for_norm(index), ndims(x))...
                     ); 
                 
% get full complex spectra
X = fft(x, [], ndims(x)) ./ N .* 2; 

% Prepare vectors with same direction as X, and unit length (make sure we
% don't divide by 0)
mX_whole_spect = abs(X);
mX_whole_spect(mX_whole_spect==0) = 1; 
X_norm_vecs = X ./ mX_whole_spect; 

% Scale them with the estimated 1/f magnitudes 
X_ap_vecs = X_norm_vecs .* ap_whole_spect; 

% prepare X normalized 
X_norm = X; 

% set X bins where the signal is BELOW noise amplitude to have amplitude of
% zero and don't touch them again
mask_dont_touch = abs(X) - abs(X_ap_vecs) < 0; 
X_norm(mask_dont_touch) = 0; 

% now wherever the magnitude of the observed spectrum is ABOVE the
% estimated 1/f noise magnitude, shorten the vector at that frequency bin
% by adding a vector with opposite direction, and magnitude equal to the
% estimated 1/f noise magnitude at that frequency. 
X_norm(~mask_dont_touch) = X(~mask_dont_touch) - X_ap_vecs(~mask_dont_touch); 

% If we know which frequency bins the signal is going to project to, we can
% simply ONLY RETAIN SIGNAL FREQUENCIES and set the complex numbers at all
% other frequency bins to zero. 
freq_to_keep_idx = [freq_to_ignore_idx; N - freq_to_ignore_idx + 2]; 

X_norm_frex_only = zeros(size(X_norm)); 
X_norm_frex_only(freq_to_keep_idx) = X_norm(freq_to_keep_idx); 


%%
    
[acf_clean, lags, ~, mX_clean] = get_acf(x_clean, fs); 

% get the autocorrelation frunciton from the normalized spectrum 
acf_subtr = real(ifft(X_norm .* conj(X_norm))); 
acf_subtr_frex_only = real(ifft(X_norm_frex_only .* conj(X_norm_frex_only))); 

% get raw autocorrelation (without 1/f normalization)
X = fft(x) / N * 2; 
acf_raw = (real(ifft(X .* conj(X)))); 

% only output acf lags up to N/2 
acf_raw = acf_raw(1:hN); 
acf_subtr = acf_subtr(1:hN); 
acf_subtr_frex_only = acf_subtr_frex_only(1:hN); 

% set acf at lag 0 to the value of lag 1 (it's just variance anyway)
if lags(1) == 0
    acf_raw(1) = acf_raw(2); 
    acf_subtr(1) = acf_subtr(2); 
    acf_subtr_frex_only(1) = acf_subtr_frex_only(2); 
end

% set DC to 0 for output (it's just the mean)
mX(1) = 0; 
mX_to_fit(1) = 0;

% make sure frequencies are a row vector
freq = ensure_row(freq);


% [mod(phase(X(1:30)), 2*pi); mod(phase(X_norm(1:30)), 2*pi); abs(X(1:30)); abs(X_norm(1:30))]'
% 
% figure
% stem(mod(phase(X(1:300)), 2*pi), 'o')
% hold on 
% stem(mod(phase(X_norm(1:300)), 2*pi), 'o')

%%

max_freq = 9;

nyq = fs/2; 
freq = [0 : hN-1]' / N * fs; 
freq_to_ignore = [response_f0 : response_f0 : nyq]'; 
freq_to_ignore_idx = dsearchn(freq, freq_to_ignore); 


f = figure('color', 'white', 'Position', [522 668 744 165]);

pnl = panel(f);

pnl.pack('v', [50, 50]);
pnl(1).pack('h', 4);
pnl(2).pack('h', 4);

for i=1:4
    pnl(2, i).pack('h', [95, 05]);
end

pnl.de.margin = [4, 3, 1, 6];
pnl.margin = [5, 2, 5, 5]; 


% ===========================================================================
% magnitude spectrum
% ===========================================================================

% clean
ax = pnl(1, 1).select(); 
plot_fft(freq, mX_clean, ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];


% raw
ax = pnl(1, 2).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];


% subtr
ax = pnl(1, 3).select(); 
plot_fft(freq, abs(X_norm(1:length(freq))), ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];

% subtr only harmonics
ax = pnl(1, 4).select(); 
plot_fft(freq, abs(X_norm_frex_only(1:length(freq))), ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];


% ===========================================================================
% acf
% ===========================================================================

% get ylimits for ACF
% ylims_acf = [...
%     min([acf_raw, acf_subtr, acf_subtr_frex_only]), ...
%     max([acf_raw, acf_subtr, acf_subtr_frex_only])...
%     ]; 

ylims_acf = [-Inf, Inf]; 

max_lag = 4.2; 

cmap = colorGradient([50, 168, 82]/255, [232, 32, 172]/255, 128); 

lags_of_interest = unique([par.lags_meter_rel, par.lags_meter_unrel]) ; 
lags_of_interest = lags_of_interest(lags_of_interest <= max_lag); 
idx = dsearchn(lags', lags_of_interest')';  


ax = pnl(2, 1, 1).select(); 

plot_acf(ax, ...
         acf_clean, ...
         lags, ...
         'min_lag', 0, ...
         'max_lag', max_lag, ...
         'lags_meter_rel', lags_of_interest, ...
         'opacity_lagz', 0.2, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;


ax = pnl(2, 1, 2).select(); 
vals = acf_clean(idx); 
imagesc(ax, 0, [1:length(vals)], vals'); 
h = pixelgrid; 
h.Children(1).Color = [1, 1, 1]; 
h.Children(1).LineWidth = 0.5; 
h.Children(2).Color = [1, 1, 1]; 
h.Children(2).LineWidth = 0.5;
colormap(ax, cmap)
ax.Visible = 'off'; 




ax = pnl(2, 2, 1).select(); 

plot_acf(ax, ...
         acf_raw, ...
         lags, ...
         'min_lag', 0, ...
         'max_lag', max_lag, ...
         'lags_meter_rel', lags_of_interest, ...
         'opacity_lagz', 0.2, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;


ax = pnl(2, 2, 2).select(); 
vals = acf_raw(idx); 
imagesc(ax, 0, [1:length(vals)], vals'); 
h = pixelgrid; 
h.Children(1).Color = [1, 1, 1]; 
h.Children(1).LineWidth = 0.5; 
h.Children(2).Color = [1, 1, 1]; 
h.Children(2).LineWidth = 0.5;
colormap(ax, cmap)
ax.Visible = 'off'; 



ax = pnl(2, 3, 1).select(); 

plot_acf(ax, ...
         acf_subtr, ...
         lags, ...
         'min_lag', 0, ...
         'max_lag', max_lag, ...
         'lags_meter_rel', lags_of_interest, ...
         'opacity_lagz', 0.2, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;

ax = pnl(2, 3, 2).select(); 
vals = acf_subtr(idx); 
imagesc(ax, 0, [1:length(vals)], vals'); 
h = pixelgrid; 
h.Children(1).Color = [1, 1, 1]; 
h.Children(1).LineWidth = 0.5; 
h.Children(2).Color = [1, 1, 1]; 
h.Children(2).LineWidth = 0.5;
colormap(ax, cmap)
ax.Visible = 'off'; 




ax = pnl(2, 4, 1).select(); 

plot_acf(ax, ...
         acf_subtr_frex_only, ...
         lags, ...
         'min_lag', 0, ...
         'max_lag', max_lag, ...
         'lags_meter_rel', lags_of_interest, ...
         'opacity_lagz', 0.2, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;


ax = pnl(2, 4, 2).select(); 
vals = acf_subtr_frex_only(idx); 
imagesc(ax, 0, [1:length(vals)], vals'); 
h = pixelgrid; 
h.Children(1).Color = [1, 1, 1]; 
h.Children(1).LineWidth = 0.5; 
h.Children(2).Color = [1, 1, 1]; 
h.Children(2).LineWidth = 0.5;
colormap(ax, cmap)
ax.Visible = 'off'; 





fname = '/datadisk/projects_backed_up/autocorrelation/figures/general/explain_ap_fit/explain_ap_fit_simple'; 

print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', f);  




%%

f = figure('color', 'white', 'Position', [587 629 400 180]);

pnl = panel(f);

pnl.pack('v', [20, 40, 40]);

pnl.de.margin = [1, 3, 1, 1];
pnl.margin = [5, 2, 5, 5]; 



ax = pnl(1).select(); 

plot(ax, t, x_clean, 'color', 'k', 'linew', 2);     
ax.XLim = [0, 3*2.4]; 
ax.Visible = 'off'; 


ax = pnl(2).select(); 

plot(ax, t, x, 'color', 'k', 'linew', 2);     
ax.XLim = [0, 3*2.4]; 
ax.Visible = 'off'; 



ax = pnl(3).select(); 

noise = get_colored_noise(length(x_clean), fs, noise_exponent); 

plot(ax, t, noise, 'color', 'k', 'linew', 2);     
ax.XLim = [0, 3*2.4]; 
ax.Visible = 'off'; 


fname = '/datadisk/projects_backed_up/autocorrelation/figures/general/explain_ap_fit/explain_ap_fit_simple_timeDomain'; 

print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', f);  





%%
% 
% 
% f = figure('color', 'white', 'Position', [587 627 151 130]);
% 
% pnl = panel(f);
% 
% pnl.pack('v',2);
% pnl(1).pack('h', [95, 5]);
% pnl(2).pack('h', [95, 5]);
% 
% ax = pnl(1, 1).select(); 
% 
% plot_acf(ax, ...
%          acf_clean, ...
%          lags, ...
%          'min_lag', 0, ...
%          'max_lag', max_lag, ...
%          'lags_meter_rel', lags_of_interest, ...
%          'opacity_lagz', 0.2, ...
%          'linew', 1, ...
%          'prec', 1e6 ...
%      ); 
% ax.YAxis.Visible = 'on';
% ax.YTick = [0]; 
% ax.XAxisLocation = 'origin'; 
% ax.YTickLabel = []; 
% 
% 
% idx = dsearchn(lags', lags_of_interest')';  
% vals = acf_clean(idx); 
% 
% ax = pnl(1, 2).select(); 
% h = imagesc(ax, 0, [1:length(vals)], vals'); 
% h = pixelgrid;
% h.Children(1).Color = [1, 1, 1]; 
% h.Children(1).LineWidth = 0.5; 
% h.Children(2).Color = [1, 1, 1]; 
% h.Children(2).LineWidth = 0.5;
% colormap(ax, brewermap([], '-PRGn'))
% ax.Visible = 'off'; 
% 
% 
% 
% acf_estim_name = 'subtrOnlyHarm'; % raw, subtr, subtrOnlyHarm
% acf_estim = acf_subtr_frex_only; % acf_raw, acf_subtr, acf_subtr_frex_only
% 
% ax = pnl(2, 1).select(); 
% 
% plot_acf(ax, ...
%          acf_estim, ... 
%          lags, ...
%          'min_lag', 0, ...
%          'max_lag', max_lag, ...
%          'lags_meter_rel', lags_of_interest, ...
%          'opacity_lagz', 0.2, ...
%          'linew', 1, ...
%          'prec', 1e6 ...
%      ); 
% ax.YAxis.Visible = 'on';
% ax.YTick = [0]; 
% ax.YTickLabel = []; 
% ax.XAxisLocation = 'origin'; 
% 
% idx = dsearchn(lags', lags_of_interest')';  
% vals = acf_estim(idx); 
% 
% ax = pnl(2, 2).select(); 
% imagesc(ax, 0, [1:length(vals)], vals'); 
% h = pixelgrid; 
% h.Children(1).Color = [1, 1, 1]; 
% h.Children(1).LineWidth = 0.5; 
% h.Children(2).Color = [1, 1, 1]; 
% h.Children(2).LineWidth = 0.5;
% colormap(ax, brewermap([], '-PRGn'))
% ax.Visible = 'off'; 
% 
% 
% pnl.de.margin = 2;
% pnl.margin =  2; 
% 
% 
% fname = sprintf('/datadisk/projects_backed_up/autocorrelation/figures/general/explain_ap_fit/explain_ap_fit_dist-%s-truth', ...
%     acf_estim_name); 
% 
% print(fname, '-dsvg', '-painters', f);  
% print(fname, '-dpng', '-painters', f);  
% 
% 
% 
% 
% 







