% This script shows that a periodically recurring signal will have FFT peaks at
% recurrence rate and harmonics. 

clear 

par = get_par(); 

addpath(genpath('lib'))
addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 

%% simulate

fs = 1000; 

    
grid_ioi = 0.2; 

ir = get_ir('square', fs); 
% 
% figure
% plot(ir)

%% 


f = figure('color','white', 'position', [197 692 1256 175]); 
pnl = panel(f); 
pnl.pack('v', 2); 
pnl(1).pack('h', [60, 10, 20]); 
pnl(2).pack('h', [60, 10, 20]); 
    

% make whole signal 
pat = repmat([1 0 0 0 ], 1, 32); 

[x, t] = get_s(...
            pat, ...
            grid_ioi, ...
            fs, ...
            'n_cycles', 1, ...
            'ir', ir ...
            );

ax = pnl(1, 1).select(); 
plot(t, x, 'linew', 3, 'color', 'k'); 
ax.XLim = [0, 5]; 
ax.XTick = [0, 2]; 
ax.Visible = 'off'; 

mX = abs(fft(x)); 
mX(1) = 0; 
N = length(x); 
freq = [0 : N-1] / N * fs; 

ax = pnl(1, 2).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'linew', 3, ...
         'frex_meter_rel', [2.5, 5], ...
         'maxfreqlim', 5); 
ax.YAxis.Visible = 'on';      
ax.YTick = []; 


[acf, lags, ~, mX, freq] = get_acf(x, fs, 'normalize_x', true);    

ax = pnl(1, 3).select(); 
plot_acf(ax, acf, lags, 'prec', 1e12, ...
    'max_lag', 2.4, ...
    'lags_meter_rel', [0.4, 1.2, 2.0], ...
    'lags_meter_unrel', [0.8 1.6], ...
    'opacity_lagz', 0.3, ...
    'linew_lagz', 3)

ax.YAxis.Visible = 'on';      
ax.YTick = []; 



% make whole signal 
pat = repmat([1 0 1 0 ], 1, 32); 

[x, t] = get_s(...
            pat, ...
            grid_ioi, ...
            fs, ...
            'n_cycles', 1, ...
            'ir', ir ...
            );

ax = pnl(2, 1).select(); 
plot(t, x, 'linew', 3, 'color', 'k'); 
ax.XLim = [0, 5]; 
ax.XTick = [0, 2]; 
ax.Visible = 'off'; 

mX = abs(fft(x)); 
mX(1) = 0; 
N = length(x); 
freq = [0 : N-1] / N * fs; 

ax = pnl(2, 2).select(); 
plot_fft(freq, mX, ...
         'ax', ax, ...
         'linew', 3, ...
         'frex_meter_rel', [2.5, 5], ...
         'maxfreqlim', 5); 
ax.YAxis.Visible = 'on';      
ax.YTick = []; 


[acf, lags, ~, mX, freq] = get_acf(x, fs, 'normalize_x', true);    

ax = pnl(2, 3).select(); 
plot_acf(ax, acf, lags, 'prec', 1e12, ...
    'max_lag', 2.4, ...
    'lags_meter_rel', [0.4, 1.2, 2.0], ...
    'lags_meter_unrel', [0.8 1.6], ...
    'opacity_lagz', 0.3, ...
    'linew_lagz', 3)

ax.YAxis.Visible = 'on';      
ax.YTick = []; 




pnl.margin = [2, 10, 0, 5]; 
pnl(1, 1).marginright = 10; 
pnl(2, 1).marginright = 10; 

pnl(1, 3).marginleft = 10; 
pnl(2, 3).marginleft = 10; 

pnl(1, 3).marginright = 0; 
pnl(2, 3).marginright = 0; 

pnl.marginright = 1; 


%%

fpath = '/datadisk/projects_backed_up/autocorrelation/figures/general/explain_beat_subdiv'; 
if ~isdir(fpath)
    mkdir(fpath)
end

fname = 'explain_beat_subdiv'; 

print(fullfile(fpath, fname), '-dsvg', '-painters', f);  
print(fullfile(fpath, fname), '-dpng', '-painters', f);  

