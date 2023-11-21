% This script shows that a periodically recurring signal will have FFT peaks at
% recurrence rate and harmonics. 

clear 

par = get_par(); 

addpath(genpath('lib'))
addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 

%% simulate

fs = 1000; 

pat = [1 1 1 1 1 1 1 1 1 1 1] % [1 1 1 0 1 1 1 0 1 1 0 0], % [1 0 1 1 1 1 0 1 1 1 0 0]
    
grid_ioi = 0.2; 

do_noise = true; 

%% 

n = 2; 

ir = get_erp_kernel(fs, ...
    'amplitudes', [0.1, 0.1, 0.2],...
    't0s', [0, 0, 0], ...
    'taus', [0.2, 0.050, 0.3], ...
    'f0s', [1, 9, 20], ...
    'duration', 0.2 ...
    ); 
% 
% figure
% plot(ir)

f = figure('color','white', 'position', [197 567 1500 130]); 
pnl = panel(f); 
pnl.pack('h', [60, 10, 20]); 

if do_noise
    jit = 0.03; 

    x = zeros(1, round(fs * length(pat) * grid_ioi)); 
    for i=1:length(pat)

        idx = max(0, round(  ((i-1)*grid_ioi + rand*jit-jit/2)  *fs)); 

        noise = rand(1, 1*fs); 
        [b, a] = butter(4, 12/(fs/2), 'low');
        noise = filtfilt(b, a, noise); 
        idx_start = length(ir)/2; 
        noise = noise(idx_start+1 : idx_start+length(ir)); 
        noise = noise - mean(noise); 
        noise = noise - noise(1); 

        x(idx+1:idx+length(ir)) = ir + 0.1*noise; 
    end
    t = [0 : length(x)-1] / fs; 

else
    
    % make whole signal 
    [x, t] = get_s(...
                pat, ...
                grid_ioi, ...
                fs, ...
                'n_cycles', 1, ...
                'ir', ir ...
                );
            
end


ax = pnl(1).select(); 
plot(t, x, 'linew', 3, 'color', 'k'); 
ax.XLim = [0, length(t)/fs]; 
ax.XTick = [0, length(t)/fs]; 
ax.Visible = 'off'; 

mX = abs(fft(x)); 
mX(1) = 0; 
N = length(x); 
freq = [0 : N-1] / N * fs; 

ax = pnl(2).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'frex_meter_rel', [1/grid_ioi : 1/grid_ioi : fs/2], ...
         'maxfreqlim', 50); 
     
ax.YLim = [0, 15]; 


[acf, lags, ~, mX, freq] = get_acf(x, fs, 'normalize_x', true);    

ax = pnl(3).select(); 
plot_acf(ax, acf, lags, 'prec', 1e12, ...
    'lags_meter_rel', [grid_ioi : grid_ioi : t(end)/2], ...
    'opacity_lagz', 0.3)

ax.YLim = [-0.0012, 0.003]; 


pnl.margin = [2, 10, 0, 5]; 
pnl(1).marginright = 10; 
pnl(3).marginleft = 10; 
pnl(3).marginright = 0; 

%%

if do_noise
    str_append = '_noise'; 
else
    str_append = ''; 
end

fname = sprintf('/datadisk/projects_backed_up/autocorrelation/figures/general/explain_fft_periodic_recurrence%s', str_append); 

print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', f);  

close(f)
