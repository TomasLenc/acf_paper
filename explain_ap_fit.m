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

f0_to_ignore = 1 / (length(pat) * grid_ioi);

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
freq_to_ignore = [f0_to_ignore : f0_to_ignore : nyq]'; 
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

% get full complex spectra
X = fft(x, [], ndims(x)) ./ N .* 2; 

% If the first frequency bin is zero, the value of estimated aperiodic
% value will be Inf. As we're going to divide by AP bin by bin, let's set
% this value to 1 (instead of Inf). This will keep the DC magnitude
% untouched. 
ap_for_norm = ap_linear; 
if freq(1) == 0
    ap_for_norm(1) = 1; 
end

% mirror the aperiodic compoent so we also have it for negative frequencies 
if mod(N, 2) == 0
    % even number of frequency bins => we have bin at pi/2 => make sure you
    % don't copy this bin twice! 
    index = [2 : length(ap_for_norm) - 1];      
    ap_whole_spect = cat(...
                         ndims(x), ...
                         ap_for_norm, ...
                         flip(ap_for_norm(index), ndims(x))...
                         ); 
else
    % odd number of frequency bins => no bin at pi/2 => just skip DC and mirror
    index = [2 :  length(ap_for_norm)];      
    ap_whole_spect = cat(...
                         ndims(x), ...
                         ap_for_norm, ...
                         flip(ap_for_norm(index), ndims(x))...
                         ); 

end

% Divide each frequency bin by the estimated 1/f magnitude at that
% frequency -> if the spectrum precisely follows the 1/f component, all
% magnitudes shold be normalized to 0. 
X_norm = X ./ ap_whole_spect; 
    
% Convert 1/f-normalized spectrum back to time domain. 
x_norm = real(ifft(X_norm)); 

% get the autocorrelation frunciton from the normalized spectrum 
acf_subtr = real(ifft(X_norm .* conj(X_norm))); 

% get raw autocorrelation (without 1/f normalization)
X = fft(x) / N * 2; 
acf_raw = (real(ifft(X .* conj(X)))); 

% only output acf lags up to N/2 
acf_raw = acf_raw(1:hN); 
acf_subtr = acf_subtr(1:hN); 

% set acf at lag 0 to the value of lag 1 (it's just variance anyway)
if lags(1) == 0
    acf_raw(1) = acf_raw(2); 
    acf_subtr(1) = acf_subtr(2); 
end

% set DC to 0 for output (it's just the mean)
mX(1) = 0; 
mX_to_fit(1) = 0;

% make sure frequencies are a row vector
freq = ensure_row(freq);


%%

max_freq = 9;


f = figure('color', 'white', 'Position', [659 92 1121 779]);

pnl = panel(f);

pnl.pack('v', [10, 20, 40, 30]);
pnl(2).pack('h', 4);
pnl(3).pack('h', 2);
pnl(4).pack('h', 2);

% pnl.select('all');

ax = pnl(1).select();

plot_erp(x, 't', t, 'col', [0, 0, 0], 'linew', 1, 'ax', ax);

ax.XLim = [0, N/fs];
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = [0, length(pat)*grid_ioi, N/fs];


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
ax = pnl(2, 3).select(); 
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



pnl(3, 1).pack('v', [70, 30]);

freq_all = [0 : N-1] / N * fs;
freq_to_ignore = [f0_to_ignore : f0_to_ignore : fs]'; 
freq_to_ignore_idx = dsearchn(ensure_col(freq_all), freq_to_ignore); 

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
plot(ax, freq_all, unwrap(phase(X)), 'linew', 0.7, 'color', [0.4, 0.4, 0.4]); 
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'on';
ax.XTick = fs * [0, 1/4, 1/2, 3/4, 1];



pnl(3, 2).pack('v', [70, 30]);

ax = pnl(3, 2, 1).select(); 
plot_fft(freq_all, abs(X_norm), ...
         'frex_meter_rel', freq_to_ignore, ...
         'ax', ax, ...
         'linew', 0.2); 
ax.XLim = [0, Inf];
ax.YAxis.Visible = 'off';

ax = pnl(3, 2, 2).select(); 
plot(ax, freq_all, unwrap(phase(X_norm)), 'linew', 0.7, 'color', [0.4, 0.4, 0.4]); 
ax.XAxis.Visible = 'on';
ax.XTick = fs * [0, 1/4, 1/2, 3/4, 1];





ax = pnl(4, 1).select(); 

plot_acf(ax, ...
         acf_raw, ...
         lags, ...
         'min_lag', 0, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';




ax = pnl(4, 2).select(); 

plot_acf(ax, ...
         acf_subtr, ...
         lags, ...
         'min_lag', 0, ...
         'linew', 1, ...
         'prec', 1e6 ...
     ); 
ax.YAxis.Visible = 'off';



pnl.de.margin = [10, 3, 5, 8];

pnl(4).margintop = 15;


%%

save_fig(f, 'explain_ap_fit');  


%%



f = figure('color', 'white', 'Position', [887 390 668 283]);

pnl = panel(f);

pnl.pack('h', 2);

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


save_fig(f, 'explain_ap_fit_complex_spectra');  















