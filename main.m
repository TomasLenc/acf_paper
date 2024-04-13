clear 

%%
[~, hostname] = system('hostname');
hostname = deblank(hostname);

par = get_par(); 

%% noies

n_noise_samples = 500;

% genarate noies
if strcmp(par.noise_type, 'eeg')

    noise_all_samples = prepare_eeg_noise(n_noise_samples, par.trial_dur);    

elseif strcmp(par.noise_type, 'fractal')

    noise_all_samples = get_colored_noise2(...
        [n_noise_samples, round(par.trial_dur*par.fs)], ...
        par.fs, par.noise_exponent); 

else
    
    error('noise type "%s" not implemented', par.noise_type);

end


%%

% % this one is independent on the frex/lags selection 
% par = get_par(); 
% 
% par.ir_type = 'erp2'; 
% par.ir = get_ir(par.ir_type, par.fs); 
% 
% par.n_rep = 50; 
% 
% par.data_path = fullfile(par.data_path); 
% 
% main_snr_vs_nlags(par,...
%     'prepared_noise', noise_all_samples); 


%% 

% sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% % frequencies of interst
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% % lags of interest 
% par.max_lag = par.trial_dur / 2; 
% 
% par.lag_base_incl_meter_rel = [0.8]; 
% par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]
% 
% par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
% par.lag_base_excl_meter_unrel = [0.8]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
%         
% run_mains


%% 
% 
% % This one is only to double check that the better robustnesss of ACF to noise
% % is not simply because we're taking more values (lags of interest). So let's
% % equalize the number of frex and lags of interest and re-run only the relevant
% % scripts. 
% 
% sel_name = 'maxlag-11lags_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% % frequencies of interst
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% % lags of interest 
% par.max_lag = par.trial_dur / 2; 
% 
% par.lag_base_incl_meter_rel = [0.8]; 
% par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]
% 
% par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
% par.lag_base_excl_meter_unrel = [0.4]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
%        
% % take only firsst 11 lags
% lags = [par.lags_meter_rel, par.lags_meter_unrel]; 
% 
% mask_meter_rel = zeros(1, length(lags), 'logical'); 
% mask_meter_rel(1 : length(par.lags_meter_rel)) = 1; 
% 
% [l, idx] = sort(lags); 
% idx = idx(1 : length(par.frex)); 
% 
% par.lags_meter_rel = lags( idx(ismember(idx, find(mask_meter_rel))) ); 
% par.lags_meter_unrel = lags( idx(ismember(idx, find(~mask_meter_rel))) ); 
%         
% par.max_lag = max([par.lags_meter_rel, par.lags_meter_unrel]) * 1.1; 
% 
% run_mains


%% 

% sel_name = 'maxlag-4.8_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% % frequencies of interst
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% % lags of interest 
% par.max_lag = 4.8; 
% 
% par.lag_base_incl_meter_rel = [0.8]; 
% par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; 
% 
% par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; 
% par.lag_base_excl_meter_unrel = [0.4]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
       
% % also run the script that direcrly compares sensitivity to noise between
% % taking the whole trial and when taking averaged erp chunks
% main_noiseEffectDist_erpVsTrial(par, 'prepared_noise', noise_all_samples);

% % alternative script to analyse lowhigh EEG data after chunking and averaging
% % erp segments within each trial - this may give higher SNR than simply taking
% % the full trial? Can compare in R. 
% main_lowhigh_eeg_chunk(par); 

% run_mains



%%

% 
% par = get_par(); 
% 
% % frequencies of interst
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% % lags of interest 
% par.max_lag = 2.4; 
% 
% par.lags_meter_rel = [0.8, 1.6]; 
% par.lags_meter_unrel = [0.6, 1.0, 1.4, 1.8]; 
% 
% 
% sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-false_keepBand-false'; 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% par.only_use_f0_harmonics = false; 
% par.ap_band_around_harmonics = [1, 1]; 
% 
% run_mains
% 
% 
% sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-false'; 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% par.only_use_f0_harmonics = true; 
% par.ap_band_around_harmonics = [1, 1]; 
% 
% run_mains
% 
% 
% sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-true'; 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% par.only_use_f0_harmonics = true; 
% par.ap_band_around_harmonics = [2, 5]; 
% 
% run_mains
% 


%% 


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




%% 

% sel_name = 'maxlag-halfTrial_meterRel-0.4_meterUnrel-0.6_1.0_1.4'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% % frequencies of interst
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% % lags of interest 
% par.max_lag = par.trial_dur / 2; 
% 
% par.lag_base_incl_meter_rel = [0.4]; 
% par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]
% 
% par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
% par.lag_base_excl_meter_unrel = [0.4]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
%         
% run_mains


%% 
% 
% sel_name = 'maxlag-halfTrial_meterUnrel-0.2'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% par.max_lag = par.trial_dur / 2; 
% 
% par.lag_base_incl_meter_rel = [0.8]; 
% par.lag_base_excl_meter_rel = [2.4]; % [0.6, 1.0, 1.4]   [2.4]
% 
% par.lag_base_incl_meter_unrel = [0.2]; % [0.6, 1.0, 1.4]   [0.2]
% par.lag_base_excl_meter_unrel = [0.8]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
%         
% run_mains
% 
% %%
% 
% sel_name = 'maxlag-2.4_meterUnrel-0.2'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
% 
% par.max_lag = 2.4; 
% 
% par.lag_base_incl_meter_rel = [0.8]; 
% par.lag_base_excl_meter_rel = [2.4]; % [0.6, 1.0, 1.4]   [2.4]
% 
% par.lag_base_incl_meter_unrel = [0.2]; % [0.6, 1.0, 1.4]   [0.2]
% par.lag_base_excl_meter_unrel = [0.8]; 
% 
% [par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
%             par.max_lag, ...
%             par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
%             par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
%             );
%         
% run_mains
% 
% 


% 
% %% 
% 
% sel_name = 'maxfreq-5_excl5'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
%         
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% 
% run_mains
% 
% 
% %% 
% 
% sel_name = 'maxfreq-30_excl5'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
%         
% par.max_freq = 30; 
% par.max_freq_plot = 30.1; 
% par.f0_to_excl = 5; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq, ...
%                                                 'f0_to_excl', par.f0_to_excl);
% 
% 
% run_mains

% 
% %% 
% 
% sel_name = 'maxfreq-5_excl0.416,0.833'; 
% 
% par = get_par(); 
% 
% par.data_path = fullfile(par.data_path, sel_name); 
% par.fig_path = par.data_path; 
% mkdir(par.data_path); 
%         
% par.max_freq = 5; 
% par.max_freq_plot = 5.1; 
% [par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
%                                                 par.max_freq);
% 
% par.freq_meter_rel(par.freq_meter_rel < 1) = []; 
% par.freq_meter_unrel(par.freq_meter_unrel < 1) = []; 
%                                             
% run_mains
% 



