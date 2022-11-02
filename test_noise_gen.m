clear 
fprintf('\n');

do_fft = true; 

%% simulate

fs = 128; 
dur = 500; 
N = round(dur * fs); 

exp_1 = -2; 
exp_2 = -2; 

x_1 = get_colored_noise(N, fs, exp_1, 'f_rotation', 1)'; 

x_2 = get_colored_noise(N, fs, exp_1, 'f_rotation', 1)'; 

% test if the multidimensional generator works
x_2 = get_colored_noise2([2, 3, N], fs, exp_1); 
x_2 = squeeze(x_2(1, 1, :))'; 


%% analyse spectra

% hN = floor(N/2)+1; 
% freq = [0:hN-1]'/N*fs; 
%
% mX_1 = abs(fft(x_1)) / hN; 
% mX_1 = mX_1(1:hN); 
% pow_1 = mX_1 .^ 2; 
% 
% mX_2 = abs(fft(x_2))/ hN; 
% mX_2 = mX_2(1:hN); 
% pow_2 = mX_2 .^ 2; 

[pow_1, freq_norm] = pwelch(x_1); 
[pow_2, freq_norm] = pwelch(x_2); 
freq = freq_norm / (2*pi) * fs; 
mX_1 = sqrt(pow_1); 
mX_2 = sqrt(pow_2); 

log_pow_1 = log(pow_1(2:end)); 
log_pow_2 = log(pow_2(2:end)); 
log_freq = log(freq(2:end)); 

log_freq_lim_min = -4; 
log_freq_lim_max = 4;  
log_freq_lim_min_idx = dsearchn(log_freq, log_freq_lim_min); 
log_freq_lim_max_idx = dsearchn(log_freq, log_freq_lim_max); 


% estimate power law from log-log
% -------------------------------

% X = [ones(size(log_freq)), log_freq]; 
% betas = (X'*X) \ (X'*log_pow_1); 
betas_1 = polyfit(log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx), ...
                log_pow_1(log_freq_lim_min_idx : log_freq_lim_max_idx), 1); 
betas_1 = flip(betas_1); 

betas_2 = polyfit(log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx), ...
                log_pow_2(log_freq_lim_min_idx : log_freq_lim_max_idx), 1); 
betas_2 = flip(betas_2); 


%% plot spectra

h = []; 

f = figure('color', 'white', 'position', [476 222 560 697]); 
pnl = panel(f); 
pnl.pack('v', 2); 

% magnitude
% ---------
pnl(1).select(); 
plot(freq, mX_1, 'linew', 2)
hold on
plot(freq, mX_2, 'linew', 2)
leg = legend({num2str(exp_1), num2str(exp_2)}); 
title(leg, 'exponent')
xlabel('frequency (Hz)')
ylabel('magnitude')

% log power
% ---------
pnl(2).select(); 
h(1) = plot(log_freq, log_pow_1, 'linew', 2); 
hold on
h(2) = plot(log_freq, log_pow_2, 'linew', 2);
xlabel('log frequency')
ylabel('log power')

y_pred = betas_1(1) + log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx) * betas_1(2); 
plot(log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx), y_pred, 'k', 'linew', 2)
plot(log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx), y_pred, 'k', 'linew', 2)

y_pred = betas_2(1) + log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx) * betas_2(2); 
plot(log_freq(log_freq_lim_min_idx : log_freq_lim_max_idx), y_pred, 'k', 'linew', 2)

leg = legend(h, {sprintf('%.3f', betas_1(2)), sprintf('%.3f', betas_2(2))}); 
title(leg, 'exponent estim')

pnl.fontsize = 12; 
pnl.margin = [20, 15, 5, 5]; 


%% 

fit_aperiodic(freq(2:end-1), log10(pow_1(2:end-1)))

fit_aperiodic(freq(2:end-1), log10(pow_2(2:end-1)))
