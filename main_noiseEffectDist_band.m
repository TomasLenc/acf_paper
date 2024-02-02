function main_noiseEffectDist_peakVsBand(par, varargin)
% 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;



%% find patterns 

par.all_pats = ...
    [1     1     0     1     1     1     0     1     0     1     0     0
     1     1     1     1     0     1     0     1     0     1     0     0
     1     1     1     0     1     0     1     1     0     1     0     0
     1     1     1     1     0     1     1     0     1     0     0     0
     1     1     1     0     1     1     0     1     0     1     0     0];

 
%% pararmeters

par.snrs = logspace(log10(0.01), log10(8), 10); 

f_step = 1 / par.trial_dur; 

f_pat = 1 / (size(par.all_pats, 2) * par.grid_ioi); 

bands = linspace(f_pat / 5, f_pat / 2, 4); 

bands_N = round(bands / f_step); 

par.bandwidths = [{[1, 1]}
    num2cell(round(bsxfun(@times, bands_N', [1/2, 1])), 2)
    ]; 


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
    'pat', 'snr', 'sample', 'bandwidth', 'r_raw', 'r_subtr', 'z_snr' ...
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

        feat_acf_orig = get_acf_features(...
                             acf_clean, lags, ...
                             par.lags_meter_rel, par.lags_meter_unrel);         

        % without aperiodic subtraction    
        [acf_raw] = get_acf(x, par.fs);    

        feat_acf_raw = get_acf_features(...
                             acf_raw, lags, ...
                             par.lags_meter_rel, par.lags_meter_unrel); 
        
        % empirical zSNR
        z_snr = get_z_snr(mX, freq, par.frex, ...
                               par.noise_bins_snr(1), ...
                               par.noise_bins_snr(2)); 

        for i_band=1:length(par.bandwidths)
            
            % get acf with noise correction using the target bandwidth
            [acf_subtracted] = get_acf(x, par.fs, ...
                   'rm_ap', true, ...
                   'ap_fit_method', par.ap_fit_method, ...
                   'f0_to_ignore', par.f0_to_ignore, ...
                   'only_use_f0_harmonics', true, ...
                   'keep_band_around_f0_harmonics', par.bandwidths{i_band} ...
                   );

            % get features
            feat_acf_subtr = get_acf_features(...
                                 acf_subtracted, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel); 


            % get measures of disctance from ground truth
            r_raw = corr(feat_acf_orig.vals', feat_acf_raw.vals'); 
            r_subtr = corr(feat_acf_orig.vals', feat_acf_subtr.vals'); 

            band_Hz = (par.bandwidths{i_band}(2) -1) * f_step; 

            rows = [...
                repmat({i_pat}, n_rep, 1), ...
                repmat({par.snrs(i_snr)}, n_rep, 1), ...
                num2cell([1:n_rep]'), ...
                repmat({band_Hz}, n_rep, 1), ...
                num2cell(r_raw'), ...
                num2cell(r_subtr'), ...
                num2cell(z_snr)...
                ];

            tbl = [tbl; rows]; 

        end
    end
end

%%

% save table
fname = sprintf('irType-%s_apFitMethod-%s_noiseEffectDistBand', ...
                par.ir_type, ...
                par.ap_fit_method); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


