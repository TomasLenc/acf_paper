function main_ir_rand(par)
% Simulates the effect of impulse response shape on different ways to
% calculate pominence of ACF at beat-related lags. 

n_rep = 100; 

%%

col_names = {
    'sample', 'zscore', 'ratio', 'contrast', 'r' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 


%%

for i_rep=1:n_rep
    
    scale = 10 * rand(); 
    offset = 100 * (rand-0.5); 
    
    ir = randn(1, floor(par.fs * par.grid_ioi)); 
    
    % make whole signal 
    [x, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );
               
    % shift and scale the whole signal 
    x = offset + scale * x; 

    % get acf withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, par.fs);   
                                                                          
    % get features
    feat_acf = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    
    
    % get mean Pearson correlation at beat-related lags                        
    r = 0;                         
    for i_lag=1:length(par.lags_meter_rel)        
        idx = round(par.lags_meter_rel(i_lag) * par.fs) + 1; 
        r = r + correlation(x, circshift(x, -idx));                         
    end    
    r = r / length(par.lags_meter_rel); 
    
    % add to table
    new_row = [...
        {i_rep}, ...
        {feat_acf.z_meter_rel}, ...
        {feat_acf.ratio_meter_rel}, ...
        {feat_acf.contrast_meter_rel}, ...
        {r} ...
        ]; 
    
    tbl = [tbl; new_row]; 
    
end

% save DATA 
fname = sprintf('ir-rand_ir');
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 



