% Try autocorrelation of different patterns. 

par = get_par; 


% pat = randsample([0, 1], 12 * 32, true); 

pat = [1 0 0 1 0 0 1 0]; 


% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', 4, ...
                    'ir', par.ir ...
                    );
                
                
    
% clean signal
[acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

plot(freq, mX_clean)

plot(lags, acf_clean)
