

clear 

par = get_par(); 

%% noies

n_rep = 50; 

noise = prepare_eeg_noise(n_rep, par.trial_dur);    


%%

sel_name = 'maxlag-4.8_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

% frequencies of interst
par.max_freq = 5; 
par.max_freq_plot = 5.1; 
par.f0_to_excl = 5; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq, ...
                                                'f0_to_excl', par.f0_to_excl);

% lags of interest 
par.max_lag = 4.8; 

par.lag_base_incl_meter_rel = [0.8]; 
par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
par.lag_base_excl_meter_unrel = [0.4]; 

[par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
            par.max_lag, ...
            par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
            par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
            );
       

%% simulate

par.noise_type = 'eeg'; % eeg, fractal

n_cond = 5; 

par.snrs = logspace(log10(0.2), log10(2), n_cond); 


%% generate signal

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', par.ir ...
                    );


%%

col_names = {
    'snr', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'ap_offset', 'ap_exponent', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

data_to_plot = []; 

for i_cond=1:n_cond
    
    snr = par.snrs(i_cond); 

    fprintf('calculating snr %d/%d\n', i_cond, n_cond)
    
    % scale the noise to the correct SNR 
    x = add_signal_noise(repmat(x_clean, n_rep, 1), noise, snr);
    
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
    par_ap = cell(n_rep, 1); 
    
    parfor i_rep=1:n_rep
        [acf_subtracted(i_rep, :), ~, ap(i_rep, :), ~, ~, par_ap(i_rep)] = ...
                                    get_acf(x(i_rep, :), par.fs, ...
                                           'rm_ap', true, ...
                                           'ap_fit_method', par.ap_fit_method, ...
                                           'f0_to_ignore', par.f0_to_ignore, ...
                                           'ap_fit_flims', par.ap_fit_flims); 
    end
    

            
    %%
    
    x_chunked = epoch_chunks(x, par.fs, 2.4 * 4); 
    
    x_erp = squeeze(mean(x_chunked, 1)); 
    
    [acf, lags, ~, mX, freq] = get_acf(x_erp, par.fs);    
    [acf_subtr, ~, ap] = get_acf(x_erp, par.fs, 'rm_ap', true, 'ap_fit_method', 'irasa');    
    
    i = 7; 
    figure
    subplot 211
    plot(freq, mX(i,:)); 
    hold on 
    plot(freq, ap(i,:), 'r--', 'linew', 2); 
    subplot 212
    plot(lags, acf(i,:)); 
    
    get_acf(x_erp(i,:), par.fs, 'rm_ap', true, 'ap_fit_method', 'irasa', 'plot_diagnostic', true);    

    
    
    
%%    
    
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
        repmat({snr}, n_rep, 1), ...
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
    
      
    data_to_plot(i_cond).snr = snr; 
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

% save  
fname = sprintf('ir-%s_noise-%s_apFitMethod-%s_onlyHarm-%s_snr', ...
               par.ir_type, ...
               par.noise_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 

save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 


