clear 

%%
[~, hostname] = system('hostname');
hostname = deblank(hostname);

par = get_par(); 

%% prepare noise

n_noise_samples = 500;

fname_noise = fullfile(par.data_path, ...
                       sprintf('n-%d_noise.mat', n_noise_samples)); 


if ~exist(fname_noise)
    
    % genarate noise
    if strcmp(par.noise_type, 'eeg')

        noise_all_samples = prepare_eeg_noise(n_noise_samples, par.trial_dur);    

    elseif strcmp(par.noise_type, 'fractal')

        noise_all_samples = get_colored_noise2(...
            [n_noise_samples, round(par.trial_dur*par.fs)], ...
            par.fs, par.noise_exponent); 

    else

        error('noise type "%s" not implemented', par.noise_type);

    end

    % save noise to disk so we do'nt have to keep re-builting it every time
    save(fname_noise, 'noise_all_samples', 'n_noise_samples'); 

else
    
    load(fname_noise)
    
end


%%

% these scripts run independent on the frex/lags selection 
par = get_par(); 

par.ir_type = 'erp2'; 
par.ir = get_ir(par.ir_type, par.fs); 
par.n_rep = 50; 
par.data_path = fullfile(par.data_path); 

main_snr_vs_nlags(par,...
    'prepared_noise', noise_all_samples); 



%% max lag = half trial duration 


par = get_par(); 

data_path = par.data_path; 

% frequencies of interst
par.max_freq = 5; 
par.max_freq_plot = 5.1; 
par.f0_to_excl = 5; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq, ...
                                                'f0_to_excl', par.f0_to_excl);

% lags of interest 
par.max_lag = par.trial_dur / 2; 

par.lag_base_incl_meter_rel = [0.8]; 
par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
par.lag_base_excl_meter_unrel = [0.4]; 

[par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
            par.max_lag, ...
            par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
            par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
            );
        

sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-false_keepBand-false'; 
par.data_path = fullfile(data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

par.only_use_f0_harmonics = false; 
par.ap_band_around_harmonics = [1, 1]; 

run_mains


sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-false'; 
par.data_path = fullfile(data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [1, 1]; 

run_mains


sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-true'; 
par.data_path = fullfile(data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [2, 5]; 

run_mains



