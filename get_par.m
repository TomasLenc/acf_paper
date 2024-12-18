function par = get_par()

[~, hostname] = system('hostname');
hostname = deblank(hostname);

if strcmpi(hostname, 'tux')
    
    acf_tools_path = '/datadisk/projects_git_dl/acf_tools/src'; 
    rnb_tools_path = '/datadisk/projects_git_dl/rnb_tools/src'; 
    lw_path = '/datadisk/projects_backed_up/autocorrelation/lib_external/letswave6'; 
    experiment_path = '/datadisk/projects_backed_up/autocorrelation'; 
    
elseif strcmpi(hostname, 'tomo-office-desktop')
    
    acf_tools_path = '/home/tomo/Documents/acf_tools/src'; 
    rnb_tools_path = '/home/tomo/Documents/rnb_tools/src'; 
    lw_path = '/home/tomo/Documents/letswave6'; 
    experiment_path = '/DATA2/autocorrelation'; 

elseif strcmpi(hostname, 'mac-BX23-117.local')
    
    acf_tools_path = '~/projects_git/acf_tools/src'; 
    rnb_tools_path = '~/projects_git/rnb_tools/src'; 
    lw_path = '~/projects_git/letswave6'; 
    experiment_path = '~/projects_backed_up/acf'; 
    
else

    error('host %s not set up yet...', hostname);
    
end

addpath(genpath(acf_tools_path)); 
addpath(genpath(rnb_tools_path)); 
addpath(genpath(lw_path)); 
addpath(genpath('lib'));

data_path = fullfile(experiment_path, 'results'); 
eeg_path = fullfile(experiment_path, 'data'); 
resting_eeg_path = fullfile(eeg_path, 'resting_eeg', 'ds004148'); 

if ~isdir(data_path)
    mkdir(data_path)
end

%%

fs = 200; 

pat = [1 0 1 1 1 1 0 1 1 1 0 0]; %  [1 1 1 0 1 1 1 0 1 1 0 0]  [1 0 1 1 1 1 0 1 1 1 0 0]

n_cycles = 20; 

grid_ioi = 0.2; 

trial_dur = n_cycles * length(pat) * grid_ioi; 

ir_type = 'square'; 

ir = get_ir(ir_type, fs); 

noise_type = 'eeg'; 

n_rep = 50; 


%% acf parameters 

ap_fit_method = 'irasa'; 

response_f0 = 1/2.4; 

only_use_f0_harmonics = true; 

ap_fit_flims = [0.1, 9]; 

ap_band_around_harmonics = [1, 1]; 


%% lags of interest

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
min_lag = 0;
max_lag = (grid_ioi * length(pat) * n_cycles) / 2; 

lag_base_incl_meter_rel = [0.8]; 
lag_base_excl_meter_rel = [2.4]; % [0.6, 1.0, 1.4]   [2.4]

lag_base_incl_meter_unrel = [0.2]; % [0.6, 1.0, 1.4]   [0.2]
lag_base_excl_meter_unrel = [0.8]; 


if ~exist('get_lag_harmonics', 'file')
    error('cant find function get_lag_harmonics: make sure you have added acf_tools to path...')
end

% meter-related lags 
% ------------------

% Make sure there's no overlap with muiltiples of meter-unrelated lags, and 
% also the pattern repetition period. 
lags_meter_rel = get_lag_harmonics(...
                            lag_base_incl_meter_rel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_rel ...
                            ); 

% meter-unrelated lags 
% --------------------

% Make sure there's no overlap with muiltiples of meter-related lags 
% (even 0.4 seconds!), and also the pattern repetition period. 

lags_meter_unrel = get_lag_harmonics(...
                            lag_base_incl_meter_unrel, ...
                            max_lag,...
                            'lag_harm_to_exclude', lag_base_excl_meter_unrel ...
                            ); 


% make sure one more time that there's no overlap between meter-rel and -unrel !!! 
assert(~any( min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel))) < 1e-9 ))

%% frequencies of interest

max_freq = 5;
max_freq_plot = max_freq + 0.1; 

f0_to_excl = [5]; 
[freq_meter_rel, freq_meter_unrel, frex] = get_meter_freq(...
                                               max_freq, ...
                                               'f0_to_excl', f0_to_excl);

%% 

noise_bins = [2, 5]; 

% noies bins for calculating the SNR of the raw spectra
% (harmonic-snippet-zsocre method as used by Rossion)
noise_bins_snr = [3, 13]; 


%%

save_figs = true; 

fontsize = 14; 

%% 



%% return structure 

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end


