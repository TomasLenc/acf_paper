% Script that generates plots necesary to explain the zSNR subtraction
% procedure (see Figure S4 in the paper). 

clear 

par = get_par; 


%% simulate

par.noise_type = 'eeg'; % eeg, fractal

%% generate signal

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', par.ir ...
                    );


%% prepare noise

noise = prepare_eeg_noise(1, par.trial_dur); 

%%

x = add_signal_noise(x_clean, noise, 0.4);

% withuout aperiodic subtraction    
[acf, ~, ~, mX, freq] = get_acf(x, par.fs);    


f = figure('color','white', 'position',[356 188 165 160]); 
pnl = panel(f);
pnl.pack('v', [20, 80]); 
pnl(1).pack('h', [4]); 
pnl(2).pack('h', [2]); 


ax = pnl(2, 1).select(); 

plot_fft(freq, mX, ...
         'ax', ax, ...
         'frex_meter_rel', [par.freq_meter_rel, par.freq_meter_unrel], ...
         'maxfreqlim', 4.9); 
 
ax.Visible = 'off'; 

pnl.margin = 3; 

[z_snr, mean_snip, idx_snip] = get_z_snr(mX, freq, [par.freq_meter_rel, par.freq_meter_unrel], 2, 5); 


ax = pnl(1, 1).select(); 

plot_fft(idx_snip, mean_snip, ...
         'ax', ax, ...
         'frex_meter_rel', 0); 
ax.XLim = [-6, 6];
ax.Visible = 'off'; 


ax = pnl(2, 2).select(); 

mX_subtr = subtract_noise_bins(mX, 2, 5); 

plot_fft(freq, mX_subtr, ...
         'ax', ax, ...
         'frex_meter_rel', [par.freq_meter_rel, par.freq_meter_unrel], ...
         'maxfreqlim', 4.9); 
 
ax.Visible = 'off'; 



%%

fpath = fullfile(par.fig_path, 'general', 'explain_z_snr'); 
fname = 'explain_z_snr'; 

mkdir(fpath); 

print(fullfile(fpath, fname), '-dsvg', '-painters', f);  
print(fullfile(fpath, fname), '-dpng', '-painters', f);  
