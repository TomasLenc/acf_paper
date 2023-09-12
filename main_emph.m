function main_emph(par)
% simulate effect of periodic emphasis on meter zscores 

n_cond = 6; 

par.emph_levels = linspace(0, 2, n_cond); 

%%

col_names = {
    'emphasis', 'z_meter_fft', 'z_meter_acf' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

%%

for i_cond=1:n_cond
    
    emph = par.emph_levels(i_cond); 
        
    fprintf('calculating %d/%d\n', i_cond, n_cond)
    
    % make whole signal 
    [x, t] = get_s(...
                        par.pat, ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', par.ir, ...
                        'emph_period', 4, ...
                        'emph_phase', 0, ...
                        'emph_magn', emph ...
                        );
    
    % get acf withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, par.fs);    
                                                                          
    % get features
    feat_acf = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    
                                                                  
    feat_fft = get_fft_features(mX, freq, ...
                                            par.freq_meter_rel, par.freq_meter_unrel); 
                             

    new_row = [{emph}, {feat_fft.z_meter_rel}, {feat_acf.z_meter_rel}]; 
    tbl = [tbl; new_row]; 
    
    data_to_plot(i_cond).emphasis = emph; 
    data_to_plot(i_cond).x = x; 
    data_to_plot(i_cond).t = t; 
    data_to_plot(i_cond).mX = mX; 
    data_to_plot(i_cond).freq = freq; 
    data_to_plot(i_cond).acf = acf; 
    data_to_plot(i_cond).lags = lags; 

end


fname = sprintf('ir-%s_emph', ...
               par.ir_type); 

% save DATA 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 


