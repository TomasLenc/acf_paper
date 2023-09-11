clear 

%%
[~, hostname] = system('hostname');
hostname = deblank(hostname);

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%%

if strcmpi(hostname, 'tux')

    noise_all_samples = prepare_eeg_noise(500, par.trial_dur);    

elseif strcmpi(hostname, 'tomo-office-desktop')

    noise_all_samples = prepare_eeg_noise(50, par.trial_dur);    

end

n_rep_snr = 200; 
n_rep_emph_vs_noise = 100; 
n_rep_snr_vs_nlags = 500; 
n_rep_only_noise = 1000; 


%% 

sel_name = 'maxlag-halfTrial_meterUnrel-0.6_1.0_1.4'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

par.max_lag = par.trial_dur / 2; 

par.lag_base_incl_meter_rel = [0.8]; 
par.lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

par.lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
par.lag_base_excl_meter_unrel = [0.8]; 

[par.lags_meter_rel, par.lags_meter_unrel] = get_meter_lags(...
            par.max_lag, ...
            par.lag_base_incl_meter_rel, par.lag_base_excl_meter_rel, ...
            par.lag_base_incl_meter_unrel, par.lag_base_excl_meter_unrel ...
            );
        
run_scripts



%% 

sel_name = 'maxlag-halfTrial_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 

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
        
run_scripts




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
% run_scripts
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
% run_scripts
% 
% 



%% 

sel_name = 'maxfreq-5_excl5'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 
        
par.max_freq = 5; 
par.max_freq_plot = 5.1; 
par.f0_to_excl = 5; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq, ...
                                                'f0_to_excl', par.f0_to_excl);


run_scripts


%% 

sel_name = 'maxfreq-30_excl5'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 
        
par.max_freq = 30; 
par.max_freq_plot = 30.1; 
par.f0_to_excl = 5; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq, ...
                                                'f0_to_excl', par.f0_to_excl);


run_scripts


%% 

sel_name = 'maxfreq-5_excl0.416,0.833'; 

par = get_par(); 

par.data_path = fullfile(par.data_path, sel_name); 
par.fig_path = par.data_path; 
mkdir(par.data_path); 
        
par.max_freq = 5; 
par.max_freq_plot = 5.1; 
[par.freq_meter_rel, par.freq_meter_unrel, par.frex] = get_meter_freq(...
                                                par.max_freq);

par.freq_meter_rel(par.freq_meter_rel < 1) = []; 
par.freq_meter_unrel(par.freq_meter_unrel < 1) = []; 
                                            
run_scripts




