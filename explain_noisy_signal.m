
clear

load /datadisk/projects_backed_up/autocorrelation/data/maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4/irType-square_apFitMethod-irasa_onlyHarm-true_noiseEffectZscoreACFvsFFT.mat



%%

n_pat = 2; 
pat_sel = [8, 1]; % 1 is low z, 8 is high z


f = figure('color', 'white', 'Position', [522 750 339 83]);

pnl = panel(f);

pnl.pack('v', n_pat);

for i=1:n_pat
    pnl(i).pack('h', 2);
end

pnl.de.margin = [1, 3, 1, 1];
pnl.margin = [5, 2, 5, 5]; 


snr = data_to_plot(1).snr(7); 

noise = prepare_eeg_noise(1, par.trial_dur);    


for i_plt=1:n_pat
    
    
    i_pat = pat_sel(i_plt); 
    
    
    t = data_to_plot(i_pat).t; 
    
    
    
    ax = pnl(i_plt, 1).select(); 

    x_clean = data_to_plot(i_pat).x_noisy{end-1}; 
    
    plot(ax, t, x_clean, 'color', 'k', 'linew', 1);     
    ax.XLim = [0, 4.8]; 
    ax.Visible = 'off'; 
    
    
    ax = pnl(i_plt, 2).select(); 
    
    x = add_signal_noise(x_clean, noise, snr);
    
    plot(ax, t, x, 'color', 'k', 'linew', 1);     
    ax.XLim = [0, 4.8]; 
    ax.Visible = 'off'; 
    
end

%%

fname = '/datadisk/projects_backed_up/autocorrelation/figures/general/explain_ap_fit/noiseEffectZscore_examples'; 

print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', f);  