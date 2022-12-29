
clear 

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))

fs = 200; 

duty_cycles = [0.050, 0.180]; 

save_figs = false; 

%%

f = figure('color', 'white', 'pos',  [437 549 1394 297]);
pnl = panel(); 
pnl.pack('v', 2); 

for i_cond=1:length(duty_cycles)

    ir = get_square_kernel(fs, ...
        'duration', duty_cycles(i_cond), ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 

    % make whole signal 
    [x_clean, t] = get_s(...
                        [1 0 1 1 1 1 0 1 1 1 0 0], ...
                        0.2, ...
                        fs, ...
                        'n_cycles', 16, ...
                        'ir', ir ...
                        );

    idx = dsearchn(t', 0.8); 

    r = correlation(x_clean, circshift(x_clean, idx)); 
    % corrcoef(x_clean, circshift(x_clean, idx))


    ax = pnl(i_cond).select(); 
    plot(t, x_clean + 1.8, 'linew', 4)
    hold on
    plot(t, circshift(x_clean, idx), 'linew', 4)
    xlim([0, 4.8])
    title(sprintf('Pearson r = %.2f', r)); 
    
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    
    % I think Peason's r is sensitive to duty cycle because variance of the signal
    % changes...
    fprintf('signal variance = %.2f\n', var(x_clean)); 
    
end


pnl.fontsize = 24; 
pnl(2).margintop = 20; 
pnl.margin = [5, 5, 5, 15]; 



if save_figs
   fname = sprintf('figures/raw_pearson_ir'); 
   print(fname, '-dsvg', '-painters', f);  
   print(fname, '-dpng', '-painters', '-r300', f);  
end





