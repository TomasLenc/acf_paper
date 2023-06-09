function main_snr_vs_nlags(par, varargin)
% clear
% par = get_par(); 

parser = inputParser; 

addParameter(parser, 'ir_type', 'square'); % square, erp, erp2
addParameter(parser, 'prepared_noise', []); 

ir_type = parser.Results.ir_type;
noise = parser.Results.prepared_noise;


addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

% number of simulated repetitions 
n_rep = 500; % 500

snrs = logspace(log10(0.2), log10(5), 5); 


%% prepare EEG noise

if ~isempty(noise)

    n_rep = size(noise, 1); 
    
else
 
    noise = prepare_eeg_noise(n_rep, par.trial_dur);    

end

%% simulate


if strcmp(ir_type, 'square')
    ir = get_square_kernel(par.fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
elseif strcmp(ir_type, 'erp')
    ir = get_erp_kernel(par.fs,...
        'amplitudes', 1,...
        't0s', 0, ...
        'taus', 0.050, ...
        'f0s', 7, ...
        'duration', 0.2 ...
        ); 
elseif strcmp(ir_type, 'erp2')
    ir = get_erp_kernel(par.fs,...
        'amplitudes', [0.4, 0.75],...
        't0s', [0, 0], ...
        'taus', [0.2, 0.050], ...
        'f0s', [1, 7], ...
        'duration', 0.5 ...
        ); 
end

col_names = {
    'max_lag', 'max_freq', 'snr', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 


% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );

x_clean = repmat(x_clean, n_rep, 1);               


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

    for i_snr=1:length(snrs)

        snr = snrs(i_snr); 

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
        [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr] = get_acf(x, par.fs, ...
                                           'rm_ap', true, ...
                                           'f0_to_ignore', 1/2.4, ...
                                           'min_freq', 0.1, ...
                                           'max_freq', 9);  

        feat_ap = []; 
        feat_ap.offset = cellfun(@(x) x(1), par_ap);            
        feat_ap.exponent = cellfun(@(x) x(2), par_ap);            

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
            repmat({max_lag}, n_rep, 1), ...
            repmat({max_freq}, n_rep, 1), ...
            repmat({snr}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            repmat({feat_fft_orig.z_meter_rel}, n_rep, 1), ...
            repmat({feat_acf_orig.z_meter_rel}, n_rep, 1), ...
            num2cell(feat_fft.z_snr), ...
            num2cell(feat_ap.offset), ...
            num2cell(feat_ap.exponent) ...
            ];


        tbl = [tbl; rows]; 

    end

end


% save table
fname = sprintf('irType-%s_nrep-%d_snrVsNlagsNfrex', ir_type, n_rep); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par', 'snrs', ...
    'lags_meter_rel', 'lags_meter_unrel', 'freq_meter_rel', 'freq_meter_unrel'); 



