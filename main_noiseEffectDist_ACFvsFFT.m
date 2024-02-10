function main_noiseEffectDist_ACFvsFFT(par, varargin)

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;


%% simulate

n_rep = par.n_rep; 

par.snrs = logspace(log10(0.01), log10(8), 10); 



%% find patterns 

% n_pats = 5; 
% 
% n_events = 12; 
% n_sounds = 7; 
% max_group_size = 4; 
% 
% all_good_pats = find_all_patterns(n_events, n_sounds, max_group_size); 
% 
% idx = randsample(size(all_good_pats, 1), n_pats); 
% 
% par.all_pats = all_good_pats(idx, :); 

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
    'r_fft_raw', 'r_acf_raw', 'r_fft_subtr', 'r_acf_subtr', ...
    'l2_fft_raw', 'l2_acf_raw', 'l2_fft_subtr', 'l2_acf_subtr', ...
    'z_snr' ...
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
        [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

        % withuout aperiodic subtraction    
        [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    

        mX_subtracted = subtract_noise_bins(mX, ...
                                par.noise_bins(1), par.noise_bins(2)); 

        % with aperiodic subtraction    
        [acf_subtracted] = get_acf(x, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'f0_to_ignore', par.f0_to_ignore, ...
                       'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                       'keep_band_around_f0_harmonics', par.ap_band_around_harmonics, ...
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


        % get features for the clean spectra
        feat_fft_orig = get_fft_features(mX_clean, freq,...
                                par.freq_meter_rel, par.freq_meter_unrel); 

        % get features for the raw spectra                                    
        tmp = get_fft_features(mX, freq, ...
                        par.freq_meter_rel, par.freq_meter_unrel); 
        feat_fft = []; 
        feat_fft.z_meter_rel = tmp.z_meter_rel; 
        feat_fft.vals = tmp.vals; 
        feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 

        % get features for the 1/f-subtracted spectra                                    
        feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                par.freq_meter_rel, par.freq_meter_unrel);

                            

        % get measures of disctance from ground truth
        pearson_acf = corr(feat_acf_orig.vals', feat_acf.vals'); 
        pearson_acf_subtr = corr(feat_acf_orig.vals', feat_acf_subtracted.vals'); 
        
        pearson_fft = corr(feat_fft_orig.vals', feat_fft.vals'); 
        pearson_fft_subtr = corr(feat_fft_orig.vals', feat_fft_subtracted.vals'); 
        
        l2_acf = pdist2(feat_acf_orig.vals, feat_acf.vals); 
        l2_acf_subtr = pdist2(feat_acf_orig.vals, feat_acf_subtracted.vals); 
        
        l2_fft = pdist2(feat_fft_orig.vals, feat_fft.vals); 
        l2_fft_subtr = pdist2(feat_fft_orig.vals, feat_fft_subtracted.vals); 
        
        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({par.snrs(i_snr)}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            num2cell(pearson_fft'), ...
            num2cell(pearson_acf'), ...
            num2cell(pearson_fft_subtr'), ...
            num2cell(pearson_acf_subtr'), ...
            num2cell(l2_fft'), ...
            num2cell(l2_acf'), ...
            num2cell(l2_fft_subtr'), ...
            num2cell(l2_acf_subtr'), ...
            num2cell(feat_fft.z_snr)...
            ];

        tbl = [tbl; rows]; 

    end
end

%%

% save table
fname = sprintf('irType-%s_apFitMethod-%s_onlyHarm-%s_noiseEffectDistACFvsFFT', ...
                par.ir_type, ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics)); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


