function main_noiseEffectDist_erpVsTrial(par, varargin)

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;


%% simulate

par.erp_chunk_dur = 2.4 * 4;

n_rep = par.n_rep; 

par.snrs = logspace(log10(0.01), log10(8), 10); 


%% find patterns 

par.all_pats = ...
    [1     1     0     1     1     1     0     1     0     1     0     0
     1     1     1     1     0     1     0     1     0     1     0     0
     1     1     1     0     1     0     1     1     0     1     0     0
     1     1     1     1     0     1     1     0     1     0     0     0
     1     1     1     0     1     1     0     1     0     1     0     0];


%% prepare noise
if size(noise, 1) < n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          n_rep, size(noise, 1)); 
else
    noise = noise(1:n_rep, :); 
end

%%

col_names = {
    'pat', 'snr', 'sample', ...
    'r_trial_raw', 'r_erp_raw', ...
    'r_trial_subtr', 'r_erp_subtr', ...
    'l2_trial_raw', 'l2_erp_raw', ...
    'l2_trial_subtr', 'l2_erp_subtr', ...
    'z_snr_trial' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

for i_pat=1:size(par.all_pats, 1)
    
    % make whole signal 
    [x_clean, t] = get_s(...
                        par.all_pats(i_pat, :), ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', par.ir);

    parfor i_snr=1:length(par.snrs)

        fprintf('pat %d, snr %.1f\n', i_pat, par.snrs(i_snr)); 

        % scale the noise to the correct SNR 
        x = add_signal_noise(x_clean, noise, par.snrs(i_snr));

        
        %% whole-trial 
        
        % get acf
        % -------

        % clean signal
        [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

        % withuout aperiodic subtraction    
        [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    

        % with aperiodic subtraction    
        [acf_subtracted] = get_acf(x, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'response_f0', par.response_f0, ...
                       'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                       'ap_fit_flims', par.ap_fit_flims ...
                       );

        % get features
        % ------------

        feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel);         

        feat_acf = get_acf_features(acf, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 

        % get features for the raw spectra                                    
        feat_fft = []; 
        feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 

                            
        % get measures of disctance from ground truth
        pearson_trial = corr(feat_acf_orig.vals', feat_acf.vals'); 
        pearson_trial_subtr = corr(feat_acf_orig.vals', feat_acf_subtracted.vals'); 
                
        l2_trial = pdist2(feat_acf_orig.vals, feat_acf.vals); 
        l2_trial_subtr = pdist2(feat_acf_orig.vals, feat_acf_subtracted.vals); 
                
        

        %% average n-cycle erp         
        
        
        x_chunked = epoch_chunks(x, par.fs, par.erp_chunk_dur); 

        x_erp = squeeze(mean(x_chunked, 1)); 

        % get acf
        % -------

        % withuout aperiodic subtraction    
        [acf, ~, ~, mX, ~] = get_acf(x_erp, par.fs);    

        % with aperiodic subtraction    
        [acf_subtracted] = get_acf(x_erp, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'response_f0', par.response_f0, ...
                       'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                       'ap_fit_flims', par.ap_fit_flims ...
                       );

        % get features
        % ------------

        feat_acf = get_acf_features(acf, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 

                            
        % get measures of disctance from ground truth
        pearson_erp = corr(feat_acf_orig.vals', feat_acf.vals'); 
        pearson_erp_subtr = corr(feat_acf_orig.vals', feat_acf_subtracted.vals'); 
                
        l2_erp = pdist2(feat_acf_orig.vals, feat_acf.vals); 
        l2_erp_subtr = pdist2(feat_acf_orig.vals, feat_acf_subtracted.vals); 
                
        
        
        %% 
        
        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({par.snrs(i_snr)}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            num2cell(pearson_trial'), ...
            num2cell(pearson_erp'), ...
            num2cell(pearson_trial_subtr'), ...
            num2cell(pearson_erp_subtr'), ...
            num2cell(l2_trial'), ...
            num2cell(l2_erp'), ...
            num2cell(l2_trial_subtr'), ...
            num2cell(l2_erp_subtr'), ...
            num2cell(feat_fft.z_snr)...
            ];

        tbl = [tbl; rows]; 

    end
end

%%

% save table
fname = sprintf('irType-%s_apFitMethod-%s_onlyHarm-%s_noiseEffectDistErpVsTrial', ...
                par.ir_type, ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics)); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


