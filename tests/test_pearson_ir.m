
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))

fs = 200; 

duty_cycles = [0.050, 0.180]; 


%%

f = figure('color', 'white', 'pos',  [437 549 1394 297]);
pnl = panel();

pnl.pack('h', [80, 20]); 
pnl(1).pack('v', 2); 
pnl(2).pack('v', 2); 

for i_cond=1:length(duty_cycles)

    ir = get_square_kernel(fs, ...
        'duration', duty_cycles(i_cond), ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 

    % make whole signal 
    [x, t] = get_s(...
                        [1 0 1 1 1 1 0 1 1 1 0 0], ...
                        0.2, ...
                        fs, ...
                        'n_cycles', 16, ...
                        'ir', ir ...
                        );

    idx = dsearchn(t', 0.8); 

    r = correlation(x, circshift(x, idx)); 
    % corrcoef(x_clean, circshift(x_clean, idx))

    ax = pnl(1, i_cond).select(); 
    plot(t, x + 1.8, 'linew', 4, 'color', [246, 109, 30]/255)
    hold on
    plot(t, circshift(x, idx), 'linew', 4, 'color', [176, 92, 153]/255)
    xlim([0, 4.8])
    title(sprintf('Pearson r = %.2f', r)); 
    
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    
    % I think Peason's r is sensitive to duty cycle because variance of the signal
    % changes...
    fprintf('signal variance = %.2f\n', var(x)); 
    
    
    [acf, lags, ~, mX, freq] = get_acf(x, fs);    

    feat = get_acf_features(acf, lags, 0.8, [0.6, 1.0, 1.2]); 
    
    ax = pnl(2, i_cond).select(); 
    plot_acf(ax, acf, lags, 'prec', 1e12, ...
        'lags_meter_rel', [0.8], ...
        'lags_meter_unrel', [0.6, 1.0, 1.2], ...
        'opacity_lagz', 0.3)

    ax.XLim = [0, 2.4];
    ax.XTick = []; 
    ax.YTick = []; 
    title(sprintf('z = %.2f', feat.z_meter_rel)); 
    
end


pnl.fontsize = 24; 
pnl.margin = [5, 5, 5, 15]; 


%%

fname = sprintf('/datadisk/projects_backed_up/autocorrelation/figures/general/raw_pearson/raw_pearson_ir'); 
print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', '-r300', f);  






