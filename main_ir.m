function main_ir(par)
% Simulates the effect of impulse response shape on meter zscores. 

n_cond = 6; 

if strcmp(par.ir_type, 'square')
    par.duty_cycles = linspace(0.050, 0.180, n_cond); 
end
if strcmp(par.ir_type, 'erp')
    par.duty_cycles = linspace(10, 4, n_cond); 
end

%%

col_names = {
    'duty_cycle', 'z_meter_fft', 'z_meter_acf' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 


%%

data_to_plot = []; 

for i_cond=1:n_cond
    
    duty_cycle = par.duty_cycles(i_cond); 
        
    fprintf('calculating %d/%d\n', i_cond, n_cond);

    if strcmp(par.ir_type, 'erp')
        ir = get_erp_kernel(par.fs,...
            'amplitudes', 1,...
            't0s', 0, ...
            'taus', 0.050, ...
            'f0s', duty_cycle, ...
            'duration', 0.2 ...
            ); 
    elseif strcmp(par.ir_type, 'square')
        ir = get_square_kernel(par.fs, ...
            'duration', duty_cycle, ...
            'rampon', 0, ...
            'rampoff', 0 ...
            ); 
    else
        error('ir kind not recognized'); 
    end
    
    % make whole signal 
    [x, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );

    % get acf withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, par.fs);    
                                                                          
    % get features
    feat_acf = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    
                                                                  
    feat_fft = get_fft_features(mX, freq, ...
                                            par.freq_meter_rel, par.freq_meter_unrel); 
                             

    new_row = [{duty_cycle}, {feat_fft.z_meter_rel}, {feat_acf.z_meter_rel}]; 
    tbl = [tbl; new_row]; 
    
    data_to_plot(i_cond).duty_cycle = duty_cycle; 
    data_to_plot(i_cond).x = x; 
    data_to_plot(i_cond).t = t; 
    data_to_plot(i_cond).mX = mX; 
    data_to_plot(i_cond).freq = freq; 
    data_to_plot(i_cond).acf = acf; 
    data_to_plot(i_cond).lags = lags; 
    
end


% save DATA 
fname = sprintf('ir-%s_ir', par.ir_type);

save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 



