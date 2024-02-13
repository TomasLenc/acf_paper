addpath(genpath('/datadisk/projects_git_dl/rnb_tools/src'));
addpath(genpath('/datadisk/projects_git_dl/acf_tools/src'));
addpath(genpath('lib'))

par = get_par;

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


f = figure('color', 'white', 'Position', [318 1299 1602 779]);

pnl = panel(f);

pnl.pack('v', [10, 20, 40, 30]);
pnl(2).pack('h', 6);
pnl(3).pack('h', 3);
pnl(4).pack('h', 3);

% pnl.select('all');

% ===========================================================================
% time-domain 
% ===========================================================================

ax = pnl(1).select();

plot_erp(x, 't', t, 'col', [0, 0, 0], 'linew', 1, 'ax', ax);

ax.XLim = [0, N/fs];
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = [0, length(pat)*grid_ioi, N/fs];

% ===========================================================================
% magnitude spectrum for 1/f fitting
% ===========================================================================

% plot raw FFT
ax = pnl(2, 1).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];


% plot 1/f component
ax = pnl(2, 2).select(); 
plot_fft(freq, mX_to_fit, ...
         'ax', ax, ...
         'frex_meter_rel', freq_to_ignore, ...
         'maxfreqlim', max_freq); 
ax.YAxis.Visible = 'off';
ax.YLim = [0, max(mX(freq_to_ignore_idx))];
ax.XAxis.Visible = 'on';
ax.XTick = [0, max_freq];

hold(ax, 'on');
plot(ax, freq, ap_linear, '--', 'color', 'k', 'linew', 2);    

% ===========================================================================
% full spectrum magnitude and phase 
% ===========================================================================

freq_all = [0 : N-1] / N * fs;
freq_to_ignore = [response_f0 : response_f0 : fs-fs/length(x)]'; 
freq_to_ignore_idx = dsearchn(ensure_col(freq_all), freq_to_ignore); 

% prepate the phase for plotting
pX = unwrap(phase(X));
pX_norm = phase(X_norm); 
pX_norm_frex_only = phase(X_norm_frex_only); 
% Where the normalized spectra were zeroed out, the phase will be equal to 0.
% Replace it with the value taken from the original spectra. 
pX_norm(abs(X_norm) == 0) = pX(abs(X_norm) == 0); 
pX_norm = unwrap(pX_norm); 

pX_norm_frex_only(abs(X_norm_frex_only) == 0) = pX(abs(X_norm_frex_only) == 0); 
pX_norm_frex_only = unwrap(pX_norm_frex_only); 

% get ylimits for unwrapped phase 
ylims_ph = [...
    min([pX, pX_norm, pX_norm_frex_only]), ...
    max([pX, pX_norm, pX_norm_frex_only])...
    ]; 



pnl(3, 1).pack('v', [70, 30]);

ax = pnl(3, 1, 1).select(); 
plot_fft(freq_all, abs(X), ...
         'frex_meter_rel', freq_to_ignore, ...
         'ax', ax, ...
         'linew', 0.2); 
ax.XLim = [0, Inf];
ax.YLim = [0, prctile(ap_whole_spect, 99.95)];
ax.YAxis.Visible = 'off';

hold(ax, 'on');
plot(ax, freq_all, ap_whole_spect, '--', 'color', 'k', 'linew', 2);    
  
ax = pnl(3, 1, 2).select();
plot(ax, freq_all, pX, 'linew', 0.7, 'color', [0.4, 0.4, 0.4]); 
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = fs * [0, 1/4, 1/2, 3/4, 1];
ax.YLim = ylims_ph; 



pnl(3, 2).pack('v', [70, 30]);

ax = pnl(3, 2, 1).select(); 
plot_fft(freq_all, abs(X_norm), ...
         'frex_meter_rel', freq_to_ignore, ...
         'ax', ax, ...
         'linew', 0.2); 
ax.XLim = [0, Inf];
ax.YAxis.Visible = 'off';

ax = pnl(3, 2, 2).select(); 
plot(ax, freq_all, pX_norm, 'linew', 0.7, 'color', [0.4, 0.4, 0.4]); 
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = fs * [0, 1/4, 1/2, 3/4, 1];
ax.YLim = ylims_ph; 



pnl(3, 3).pack('v', [70, 30]);

ax = pnl(3, 3, 1).select(); 
plot_fft(freq_all, abs(X_norm_frex_only), ...
         'frex_meter_rel', freq_to_ignore, ...
         'ax', ax, ...
         'linew', 0.2); 
ax.XLim = [0, Inf];
ax.YAxis.Visible = 'off';

ax = pnl(3, 3, 2).select(); 
plot(ax, freq_all, pX_norm_frex_only, 'linew', 0.7, 'color', [0.4, 0.4, 0.4]); 
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = fs * [0, 1/4, 1/2, 3/4, 1];
ax.YLim = ylims_ph; 


% ===========================================================================
% acf
% ===========================================================================

% get ylimits for ACF
ylims_acf = [...
    min([acf_raw, acf_subtr, acf_subtr_frex_only]), ...
    max([acf_raw, acf_subtr, acf_subtr_frex_only])...
    ]; 

ax = pnl(4, 1).select(); 

plot_acf(ax, ...
         acf_raw, ...
         lags, ...
         'min_lag', 0, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;

ax = pnl(4, 2).select(); 

plot_acf(ax, ...
         acf_subtr, ...
         lags, ...
         'min_lag', 0, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;

ax = pnl(4, 3).select(); 

plot_acf(ax, ...
         acf_subtr_frex_only, ...
         lags, ...
         'min_lag', 0, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';
ax.YLim = ylims_acf;


pnl.de.margin = [10, 3, 5, 8];

pnl(4).margintop = 15;




save_fig(f, 'explain_ap_fit');  


%%


f = figure('color', 'white', 'Position', [372 1555 655 647]);

pnl = panel(f);

pnl.pack('h', 1);

% pnl.select('all');

ax = pnl(1).select();
% skip DC
for i_f=1:length(freq_all)
    c = [0.7, 0.7, 0.7];
    if ismember(i_f, freq_to_ignore_idx)
        c = [222 45 38]/255;
    end
    plot3(ax, ...
        [freq_all(i_f), freq_all(i_f)], ...
        [0, real(X(i_f))],...
        [0, imag(X(i_f))], ...
        'color', c, 'linew', 2)
    hold on
    plot3(ax, ...
        [freq_all(i_f), freq_all(i_f)], ...
        [0, real(X_ap_vecs(i_f))],...
        [0, imag(X_ap_vecs(i_f))], ...
        'color', [0, 0, 0], 'linew', 2)
end
grid on
ax.XLabel.String = 'freq';
ax.YLabel.String = 'real';
ax.ZLabel.String = 'imag';
ax.YTick = [];
ax.ZTick = [];

save_fig(f, 'explain_ap_fit_complex_spectra_ap_vectors');  


%%


f = figure('color', 'white', 'Position', [887 390 1000 283]);

pnl = panel(f);

pnl.pack('h', 3);

% pnl.select('all');

ax = pnl(1).select();
for i_f=1:length(freq_all)
    c = [0.7, 0.7, 0.7];
    if ismember(i_f, freq_to_ignore_idx)
        c = [222 45 38]/255;
    end
    plot3(ax, ...
        [freq_all(i_f), freq_all(i_f)], ...
        [0, real(X(i_f))],...
        [0, imag(X(i_f))], ...
        'color', c, 'linew', 2)
    hold on
end
grid on
ax.XLabel.String = 'freq';
ax.YLabel.String = 'real';
ax.ZLabel.String = 'imag';
ax.YTick = [];
ax.ZTick = [];


ax = pnl(2).select();
for i_f=1:length(freq_all)
    c = [0.7, 0.7, 0.7];
    if ismember(i_f, freq_to_ignore_idx)
        c = [222 45 38]/255;
    end
    plot3(ax, ...
        [freq_all(i_f), freq_all(i_f)], ...
        [0, real(X_norm(i_f))],...
        [0, imag(X_norm(i_f))], ...
        'color', c, 'linew', 2)
    hold on
end
grid on
ax.XLabel.String = 'freq';
ax.YLabel.String = 'real';
ax.ZLabel.String = 'imag';
ax.YTick = [];
ax.ZTick = [];


ax = pnl(3).select();
for i_f=1:length(freq_all)
    c = [0.7, 0.7, 0.7];
    if ismember(i_f, freq_to_ignore_idx)
        c = [222 45 38]/255;
    end
    plot3(ax, ...
        [freq_all(i_f), freq_all(i_f)], ...
        [0, real(X_norm_frex_only(i_f))],...
        [0, imag(X_norm_frex_only(i_f))], ...
        'color', c, 'linew', 2)
    hold on
end
grid on
ax.XLabel.String = 'freq';
ax.YLabel.String = 'real';
ax.ZLabel.String = 'imag';
ax.YTick = [];
ax.ZTick = [];

save_fig(f, 'explain_ap_fit_complex_spectra');  




%%

xlims_start = [0,100]; 
xlims_end = [N-100,N]; 

figure('color', 'w')
ax = subplot(2,2,1); 
plot(real(X)); 
hold on 
plot(freq_to_ignore_idx, real(X(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_start); 

ax = subplot(2,2,2); 
plot(real(X)); 
hold on 
plot(freq_to_ignore_idx, real(X(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_end); 


ax = subplot(2,2,3); 
plot(imag(X)); 
hold on
plot(freq_to_ignore_idx, imag(X(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_start); 

ax = subplot(2,2,4); 
plot(imag(X)); 
hold on
plot(freq_to_ignore_idx, imag(X(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_end); 


%%

xlims_start = [0,100]; 
xlims_end = [N-100,N]; 

figure('color', 'w')
ax = subplot(2,2,1); 
plot(real(X_ap_vecs)); 
hold on 
plot(freq_to_ignore_idx, real(X_ap_vecs(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_start); 

ax = subplot(2,2,2); 
plot(real(X_ap_vecs)); 
hold on 
plot(freq_to_ignore_idx, real(X_ap_vecs(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_end); 


ax = subplot(2,2,3); 
plot(imag(X_ap_vecs)); 
hold on
plot(freq_to_ignore_idx, imag(X_ap_vecs(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_start); 

ax = subplot(2,2,4); 
plot(imag(X_ap_vecs)); 
hold on
plot(freq_to_ignore_idx, imag(X_ap_vecs(freq_to_ignore_idx)), 'ro', 'markerfacecolor', 'r', 'MarkerSize', 2); 
box off
ax.YAxis.Visible = 'off'; 
xlim(xlims_end); 



























