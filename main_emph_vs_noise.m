function main_emph_vs_noise(par, varargin)
% Simulate the effect of noise on meter zscores, depending on the periodicity
% of the original (clean) signal. We will simulate signals based on a range of
% patterns (some very strongly periodic, some weakly periodic). On top of that,
% we will add various amounts of periodic emphais to the simulated signals for
% each pattern, resulting in a set of signals with a wide range of periodicities. 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); 

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;

%% simulate

par.emph_levels = 0; % linspace(0, 2, 4)

par.snrs = logspace(log10(0.01), log10(8), 10); 

% par.pats = {
%    % 1       1       1    
%     [1 1 1 0 1 1 1 0 1 1 0 0]
%     [1 0 1 1 1 1 0 1 1 1 0 0]
%     [1 0 1 0 1 1 0 1 1 1 0 1]
%     [1 1 0 1 1 0 0 1 1 0 1 1]
% };

par.pats = {
    [1 0 1 1 1 1 0 0 1 0 1 0]
    [1 0 1 1 1 1 0 0 1 0 0 1]
    [1 1 1 1 0 0 0 1 1 0 1 0]
    [1 1 1 0 1 1 1 0 1 0 0 0]
    [1 1 1 1 0 1 1 0 1 0 0 0]
    [1 1 1 0 1 1 0 0 1 0 1 0]
    [1 1 1 0 1 0 1 0 1 1 0 0]
    [1 1 1 0 1 0 0 1 1 0 1 0]
    [1 1 1 0 1 0 0 1 0 1 1 0]
    [1 0 1 1 1 0 0 1 1 1 0 0]
    [1 0 1 1 1 0 0 1 1 0 0 1]
    [1 1 1 0 0 1 0 1 1 0 1 0]
    [1 1 0 1 1 1 0 0 1 0 1 0]
    [1 1 0 1 1 0 0 1 1 0 1 0]
    [1 1 0 1 0 1 1 0 1 0 1 0]
    [1 0 0 1 0 1 1 1 1 0 1 0]
};


%% prepare noise

if size(noise, 1) < par.n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          par.n_rep, size(noise, 1)); 
else
    noise = noise(1:par.n_rep, :); 
end


%%

col_names = {
    'pat', 'emph', 'snr', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

for i_pat=1:length(par.pats)

    for i_emph=1:length(par.emph_levels)

        emph = par.emph_levels(i_emph); 

        % make whole signal 
        [x_clean, t] = get_s(...
                            par.pats{i_pat}, ...
                            par.grid_ioi, ...
                            par.fs, ...
                            'n_cycles', par.n_cycles, ...
                            'ir', par.ir, ...
                            'emph_period', 4, ...
                            'emph_phase', 0, ...
                            'emph_magn', emph ...
                            );

        for i_snr=1:length(par.snrs)

            snr = par.snrs(i_snr); 

            fprintf('pat %d emph %.1f, snr %.1f\n', i_pat, emph, snr); 

            % scale the noise to the correct SNR 
            x = add_signal_noise(x_clean, noise, snr);

            % get acf
            % -------

            % clean signal
            [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

            % withuout aperiodic subtraction    
            [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    

            mX_subtracted = subtract_noise_bins(mX, ...
                                    par.noise_bins(1), par.noise_bins(2)); 

            % with aperiodic subtraction    
            acf_subtracted = nan(size(acf)); 
            ap = nan(size(acf)); 
            par_ap = cell(par.n_rep, 1); 

            parfor i_rep=1:par.n_rep
                [acf_subtracted(i_rep, :), ~, ap(i_rep, :), ~, ~, ~, ~, par_ap(i_rep)] = ...
                                            get_acf(x(i_rep, :), par.fs, ...
                                                   'rm_ap', true, ...
                                                   'ap_fit_method', par.ap_fit_method, ...
                                                   'f0_to_ignore', par.f0_to_ignore, ...
                                                   'ap_fit_flims', par.ap_fit_flims, ...
                                                   'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                                                   'keep_band_around_f0_harmonics', par.ap_band_around_harmonics ...
                                                   ); 
            end

            feat_ap = []; 
            if strcmp(par.ap_fit_method, 'fooof')
                feat_ap.offset = cellfun(@(x) x(1), par_ap);            
                feat_ap.exponent = cellfun(@(x) x(2), par_ap);            
            else
                feat_ap.offset = nan(par.n_rep, 1);        
                feat_ap.exponent = nan(par.n_rep, 1);         
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
            tmp = get_fft_features(mX, freq, ...
                            par.freq_meter_rel, par.freq_meter_unrel); 
            feat_fft = []; 
            feat_fft.z_meter_rel = tmp.z_meter_rel; 

            feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                               par.noise_bins_snr(1), ...
                                               par.noise_bins_snr(2)); 

            % get features for the 1/f-subtracted spectra                                    
            feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                    par.freq_meter_rel, par.freq_meter_unrel);

            rows = [...
                repmat({i_pat}, par.n_rep, 1), ...
                repmat({emph}, par.n_rep, 1), ...
                repmat({snr}, par.n_rep, 1), ...
                num2cell([1:par.n_rep]'), ...
                num2cell(feat_fft.z_meter_rel), ...
                num2cell(feat_acf.z_meter_rel), ...
                num2cell(feat_fft_subtracted.z_meter_rel), ...
                num2cell(feat_acf_subtracted.z_meter_rel), ...
                repmat({feat_fft_orig.z_meter_rel}, par.n_rep, 1), ...
                repmat({feat_acf_orig.z_meter_rel}, par.n_rep, 1), ...
                num2cell(feat_fft.z_snr), ...
                num2cell(feat_ap.offset), ...
                num2cell(feat_ap.exponent) ...
                ];


            tbl = [tbl; rows]; 

        end
    end
end


% save table
fname = sprintf('ir-%s_noise-%s_apFitMethod-%s_onlyHarm-%s_emphVsNoise', ...
               par.ir_type,...
               par.noise_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 
           
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


