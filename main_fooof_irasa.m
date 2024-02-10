function main_fooof_irasa(par, varargin)
% Simulate signals with a range of noise. Test how much the noise distorts the
% ACF estaimates obtained with eiher fooof or irasa method to correct for the
% 1/f. The whole ground truth autocorrelation function is correlated with the
% estimated ACF using the two methods. We also test the approach where after
% subtracting the estimated 1/f component, only harmonics of the rhythmic
% pattern are taken into account when calculating the ACF. 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []);

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;

%% simulate

% number of simulated repetitions 
par.n_rep = 200; 

n_cond = 5; 

if strcmp(par.noise_type, 'eeg')
    par.snrs = logspace(log10(0.2), log10(2), n_cond); 
else
    par.snrs = logspace(log10(0.2), log10(2), n_cond); 
end

%% generate signal

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', par.ir ...
                    );

%% genearet noise

if ~isempty(noise)

    if par.n_rep ~= size(noise, 1)
        warning('requested number of repetitions smaller than provided noise...updating n_reps to %d', par.n_rep); 
        par.n_rep = size(noise, 1); 
    end
    
else
    
    if strcmp(par.noise_type, 'fractal')

        noise = get_colored_noise2([par.n_rep, length(x_clean)], par.fs, par.noise_exponent); 

    elseif strcmp(par.noise_type, 'eeg')

        noise = prepare_eeg_noise(par.n_rep, par.trial_dur);    

    else
        error('noise type "%s" not implemented', par.noise_type);
    end
    
end

%% run

col_names = {'snr', 'only_f0_harmonics', 'r_fooof', 'r_irasa'};

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


for only_use_f0_harmonics=[true, false]
    
    rows_all = cell(1, n_cond); 

    parfor i_cond=1:n_cond

        snr = par.snrs(i_cond); 

        fprintf('calculating snr %d/%d\n', i_cond, n_cond)

        % scale the noise to the correct SNR 
        x = add_signal_noise(repmat(x_clean, par.n_rep, 1), noise, snr);

        % get acf
        % -------

        % clean signal
        [acf_clean, lags] = get_acf(x_clean, par.fs);    

        % with aperiodic subtraction using FOOOF
        [acf_fooof] = get_acf(x, par.fs, ...
                               'rm_ap', true, ...
                               'ap_fit_method', 'fooof', ...
                               'f0_to_ignore', par.f0_to_ignore, ...
                               'only_use_f0_harmonics', only_use_f0_harmonics, ...
                               'keep_band_around_f0_harmonics', par.ap_band_around_harmonics, ...
                               'ap_fit_flims', par.ap_fit_flims);  

        % with aperiodic subtraction using IRASA
        [acf_irasa] = get_acf(x, par.fs, ...
                               'rm_ap', true, ...
                               'ap_fit_method', 'irasa', ...
                               'f0_to_ignore', par.f0_to_ignore, ...
                               'only_use_f0_harmonics', only_use_f0_harmonics, ...
                               'keep_band_around_f0_harmonics', par.ap_band_around_harmonics, ...
                               'verbose', 1);  

        feat_acf_clean = get_acf_features(acf_clean, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel);         

        feat_acf_fooof = get_acf_features(acf_fooof, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_irasa = get_acf_features(acf_irasa, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 
                                                      
        % get correelation between ground truth ACF and noise-subtracted one
        r_fooof = correlation(repmat(feat_acf_clean.vals, par.n_rep, 1), ...
                              feat_acf_fooof.vals); 
                          
        r_irasa = correlation(repmat(feat_acf_clean.vals, par.n_rep, 1), ...
                              feat_acf_irasa.vals); 

        rows_all{i_cond} = [...
            repmat({snr}, par.n_rep, 1), ...
            repmat({only_use_f0_harmonics}, par.n_rep, 1), ...
            num2cell(r_fooof), ...
            num2cell(r_irasa) ...
            ];               
    end
    
    % update table 
    rows = vertcat(rows_all{:});

    tbl = [tbl; rows];

end

fname = sprintf('ir-%s_noise-%s_irasaVsFooof', ...
               par.ir_type, par.noise_type); 

writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 












