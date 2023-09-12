function main_snr_vs_nlags(par, varargin)
% This script simulates how the number of frequencies/lags of interest makes
% the FFT/ACF method more or less robust to noise. 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); 

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;

par.snrs = logspace(log10(0.2), log10(5), 5); 

par.freq_meter_rel = [];
par.freq_meter_unrel = [];
par.lags_meter_rel = [];
par.lags_meter_unrel = [];


%% prepare EEG noise

if size(noise, 1) < par.n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          par.n_rep, size(noise, 1)); 
else
    noise = noise(1:par.n_rep, :); 
end

%% simulate

col_names = {
    'max_lag', 'max_freq', 'snr', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'z_snr' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', par.ir ...
                    );

x_clean = repmat(x_clean, par.n_rep, 1);               

for i_lag_sel=1:2

    min_lag = 0;

    if i_lag_sel == 1

        max_lag = 2.4; 
        max_freq = 5; 

        freq_meter_rel = 1/2.4 * [1 : 12];
        freq_meter_unrel = get_lag_harmonics(1/2.4, 5, 'lag_harm_to_exclude', 1.25);


    elseif i_lag_sel == 2

        max_lag = par.trial_dur / 2; 
        max_freq = 30; 

        freq_meter_rel = 1/2.4 * [1 : floor(max_freq/(1/2.4))];
        freq_meter_unrel = get_lag_harmonics(1/2.4, 30, 'lag_harm_to_exclude', 1.25);

    end

    lags_meter_rel = get_lag_harmonics(...
                                par.lag_base_incl_meter_rel, ...
                                max_lag,...
                                'lag_harm_to_exclude', par.lag_base_excl_meter_rel ...
                                ); 

    lags_meter_unrel = get_lag_harmonics(...
                                par.lag_base_incl_meter_unrel, ...
                                max_lag,...
                                'lag_harm_to_exclude', par.lag_base_excl_meter_unrel ...
                                ); 

    % make sure one more time that there's no overlap between meter-rel and -unrel !!!
    assert(~any( min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel))) < 1e-9 ))

    for i_snr=1:length(par.snrs)

        snr = par.snrs(i_snr); 

        fprintf('max lag %.1f, snr %.1f\n', max_lag, snr); 

        % scale the noise to the correct SNR 
        x = add_signal_noise(x_clean, noise, snr);

        % get acf
        % -------

        % clean signal
        [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean(1,:), par.fs);    

        % withuout aperiodic subtraction    
        [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    

        mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 

        % with aperiodic subtraction    
        acf_subtracted = nan(size(acf)); 
        
        parfor i_rep=1:par.n_rep
            [acf_subtracted(i_rep, :)] = ...
                                get_acf(x(i_rep, :), par.fs, ...
                                       'rm_ap', true, ...
                                       'ap_fit_method', par.ap_fit_method, ...
                                       'f0_to_ignore', par.f0_to_ignore, ...
                                       'ap_fit_flims', par.ap_fit_flims); 
        end
    
        % get features
        % ------------

        feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                     lags_meter_rel, lags_meter_unrel);         

        feat_acf = get_acf_features(acf, lags, ...
                                    lags_meter_rel, lags_meter_unrel);    

        feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                     lags_meter_rel, lags_meter_unrel); 

        % get features for the clean spectra
        feat_fft_orig = get_fft_features(mX_clean, freq,...
                                                freq_meter_rel, freq_meter_unrel); 

        % get features for the raw spectra                                    
        tmp = get_fft_features(mX, freq, freq_meter_rel, freq_meter_unrel); 
        feat_fft = []; 
        feat_fft.z_meter_rel = tmp.z_meter_rel; 

        feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 

        % get features for the 1/f-subtracted spectra                                    
        feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                               freq_meter_rel, freq_meter_unrel);

        rows = [...
            repmat({max_lag}, par.n_rep, 1), ...
            repmat({max_freq}, par.n_rep, 1), ...
            repmat({snr}, par.n_rep, 1), ...
            num2cell([1:par.n_rep]'), ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            repmat({feat_fft_orig.z_meter_rel}, par.n_rep, 1), ...
            repmat({feat_acf_orig.z_meter_rel}, par.n_rep, 1), ...
            num2cell(feat_fft.z_snr) ...
            ];

        tbl = [tbl; rows]; 

    end

    par.freq_meter_rel{i_lag_sel} = freq_meter_rel; 
    par.freq_meter_unrel{i_lag_sel} = freq_meter_unrel; 
    par.lags_meter_rel{i_lag_sel} = lags_meter_rel; 
    par.lags_meter_unrel{i_lag_sel} = lags_meter_unrel; 
    
end


fname = sprintf('ir-%s_noise-%s_apFitMethod-%s_onlyHarm-%s_snrVsNlagsNfrex', ...
               par.ir_type,...
               par.noise_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 

writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 



