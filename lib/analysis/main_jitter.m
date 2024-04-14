function main_jitter(par)
% Simluates the effect of jitter on meter zscore. 

par.pat = [1 0 0 0 1 0 0 0 1 0 0 0 ]; 

n_cond = 5; 

par.jitters = logspace(log10(0.001), log10(0.200), n_cond); 

n_rep = par.n_rep; 

%%

col_names = {
    'jitter', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'ap_offset', 'ap_exponent', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 


data_to_plot = []; 

%%

x_clean = get_s(...
                par.pat, ...
                par.grid_ioi, ...
                par.fs, ...
                'n_cycles', par.n_cycles, ...
                'ir', par.ir...
                );


%%

for i_cond=1:n_cond
    
    jit = par.jitters(i_cond); 
        
    fprintf('calculating %d/%d\n', i_cond, n_cond)
    
    % make whole signal 
    x = nan(n_rep, round(par.n_cycles * par.grid_ioi * length(par.pat) * par.fs)); 
    
    for i_rep=1:n_rep
        [x(i_rep, :), t] = get_s(...
                            par.pat, ...
                            par.grid_ioi, ...
                            par.fs, ...
                            'n_cycles', par.n_cycles, ...
                            'ir', par.ir, ...
                            'jitter', jit ...
                            );
    end
    
    % get acf
    % -------
                                       
    % clean signal
    [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

    % withuout aperiodic subtraction    
    [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    
                                   
    mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1), par.noise_bins(2)); 
    
    % with aperiodic subtraction    
    acf_subtracted = nan(size(acf)); 
    ap = nan(size(acf)); 
    
    parfor i_rep=1:n_rep
        [acf_subtracted(i_rep, :), ~, ap(i_rep, :), ~, ~, ~, ~, par_ap(i_rep)] = ...
                                    get_acf(x(i_rep, :), par.fs, ...
                                           'rm_ap', true, ...
                                           'ap_fit_method', par.ap_fit_method, ...
                                           'response_f0', par.response_f0, ...
                                           'ap_fit_flims', par.ap_fit_flims, ...
                                           'keep_band_around_f0_harmonics', par.ap_band_around_harmonics); 
    end
                                       
    feat_ap = []; 
    if strcmp(par.ap_fit_method, 'fooof')
        feat_ap.offset = cellfun(@(x) x(1), par_ap);            
        feat_ap.exponent = cellfun(@(x) x(2), par_ap);            
    else
        feat_ap.offset = nan(n_rep, 1);        
        feat_ap.exponent = nan(n_rep, 1);         
    end
    
    % get features
    % ------------
    feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

    feat_acf = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    

    feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel); 


    % get features for the clean spectra
    feat_fft_orig = get_fft_features(mX_clean, freq,...
                                            par.freq_meter_rel, par.freq_meter_unrel);
                                        
    % get features for the raw spectra                                    
    tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
    
    feat_fft = []; 
    
    feat_fft.z_meter_rel = tmp.z_meter_rel; 
                                        
    feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                       par.noise_bins_snr(1), ...
                                       par.noise_bins_snr(2)); 

    % get features for the 1/f-subtracted spectra                                    
    feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                           par.freq_meter_rel, par.freq_meter_unrel);
    
                                       
    % add features to table
    rows = [...
        repmat({jit}, n_rep, 1), ...
        num2cell([1 : n_rep]'), ...
        num2cell(feat_fft.z_meter_rel), ...
        num2cell(feat_acf.z_meter_rel), ...
        num2cell(feat_fft_subtracted.z_meter_rel), ...
        num2cell(feat_acf_subtracted.z_meter_rel), ...
        repmat({feat_fft_orig.z_meter_rel}, n_rep, 1), ...
        repmat({feat_acf_orig.z_meter_rel}, n_rep, 1), ...
        num2cell(feat_ap.offset), ...
        num2cell(feat_ap.exponent), ...
        num2cell(feat_fft.z_snr) ...
        ];

    tbl = [tbl; rows];
    
    
    data_to_plot(i_cond).jitter = jit; 
    data_to_plot(i_cond).x = x(1, :); 
    data_to_plot(i_cond).t = t; 
    data_to_plot(i_cond).mX = mX(1, :); 
    data_to_plot(i_cond).mX_subtr = mX_subtracted(1, :); 
    data_to_plot(i_cond).freq = freq; 
    data_to_plot(i_cond).acf = acf(1, :); 
    data_to_plot(i_cond).acf_subtr = acf_subtracted(1, :); 
    data_to_plot(i_cond).ap = ap(1, :); 
    data_to_plot(i_cond).lags = lags; 
    

end

%%


fname = sprintf('ir-%s_apFitMethod-%s_onlyHarm-%s_jitter', ...
               par.ir_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 

% save DATA 
save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 


