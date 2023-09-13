function main_noiseEffectDist_allVsOnlyHarm(par, varargin)
% Simulate the correlation between ACF values at lags of interest taken from
% a noise-free ground truth version of a signal, and noisy version. The idea is
% to check whether only taking the pattern harmonics when calculating the ACF
% distorts the ACF, or whether it improves the estimate. SO we will compare the
% correlation with the ground trush ACF, after estimating ACF from a noisy
% signal after (a) only subtracting the 1/f, (b) subtracting 1/f AND only
% taking harmonics of the pattern repetition rate to build the ACF. 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;


%% simulate

par.snrs = logspace(log10(0.01), log10(8), 10); 



%% find patterns 


par.all_pats = ...
    [1     1     0     1     1     1     0     1     0     1     0     0
     1     1     1     1     0     1     0     1     0     1     0     0
     1     1     1     0     1     0     1     1     0     1     0     0
     1     1     1     1     0     1     1     0     1     0     0     0
     1     1     1     0     1     1     0     1     0     1     0     0];

%% prepare noise

n_rep = par.n_rep; 

if size(noise, 1) < n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          n_rep, size(noise, 1)); 
else
    noise = noise(1:n_rep, :); 
end

%%

col_names = {
    'pat', 'snr', 'sample', 'method', 'r', 'z_snr' ...
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

        % get acf
        % -------

        % clean signal
        [acf_clean, lags] = get_acf(x_clean, par.fs);    

        % without aperiodic subtraction    
        [acf_raw] = get_acf(x, par.fs);    

        % with aperiodic subtraction    
        [acf_subtracted, ~, ~, mX, freq] = get_acf(x, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'f0_to_ignore', par.f0_to_ignore, ...
                       'only_use_f0_harmonics', false, ...
                       'ap_fit_flims', par.ap_fit_flims ...
                       );

        % with aperiodic subtraction ONLY HARMONICS
        [acf_subtracted_only_harm] = get_acf(x, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'f0_to_ignore', par.f0_to_ignore, ...
                       'only_use_f0_harmonics', true, ...
                       'ap_fit_flims', par.ap_fit_flims ...
                       );
                   
        % get features
        % ------------

        feat_acf_orig = get_acf_features(...
                                     acf_clean, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel);         

        feat_acf_raw = get_acf_features(...
                                     acf_raw, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 
                                 
        feat_acf_subtr = get_acf_features(...
                                     acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 
                                 
        feat_acf_subtr_only_harm = get_acf_features(...
                                     acf_subtracted_only_harm, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 

        z_snr = get_z_snr(mX, freq, par.frex, ...
                               par.noise_bins_snr(1), ...
                               par.noise_bins_snr(2)); 

        % get measures of disctance from ground truth
        r_raw = corr(feat_acf_orig.vals', feat_acf_raw.vals'); 
        r_subtr = corr(feat_acf_orig.vals', feat_acf_subtr.vals'); 
        r_subtr_only_harm = corr(feat_acf_orig.vals', feat_acf_subtr_only_harm.vals'); 
                
        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({par.snrs(i_snr)}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            repmat({'raw'}, n_rep, 1), ...
            num2cell(r_raw'), ...
            num2cell(z_snr)...
            ];

        tbl = [tbl; rows]; 
        
        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({par.snrs(i_snr)}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            repmat({'subtr'}, n_rep, 1), ...
            num2cell(r_subtr'), ...
            num2cell(z_snr)...
            ];

        tbl = [tbl; rows]; 

        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({par.snrs(i_snr)}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            repmat({'subtr_only_harm'}, n_rep, 1), ...
            num2cell(r_subtr_only_harm'), ...
            num2cell(z_snr)...
            ];

        tbl = [tbl; rows]; 
        
    end
end

%%

% save table
fname = sprintf('irType-%s_apFitMethod-%s_noiseEffectDistAllVsOnlyHarm', ...
                par.ir_type, ...
                par.ap_fit_method); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


