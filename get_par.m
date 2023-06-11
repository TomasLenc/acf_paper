function par = get_par()

[~, hostname] = system('hostname');
hostname = deblank(hostname);

if strcmpi(hostname, 'tux')
    
    acf_tools_path = '/datadisk/projects_git_dl/acf_tools/src'; 
    rnb_tools_path = '/datadisk/projects_git_dl/rnb_tools/src'; 
    lw_path = '/datadisk/projects_backed_up/autocorrelation/lib_external/letswave6'; 
    pica_path = '/datadisk/projects_backed_up/autocorrelation/lib_external/piCA'; 
    experiment_path = '/datadisk/projects_backed_up/autocorrelation'; 
    
elseif strcmpi(hostname, 'tomo-office-desktop')
    
    acf_tools_path = '/home/tomo/Documents/acf_tools/src'; 
    rnb_tools_path = '/home/tomo/Documents/rnb_tools/src'; 
    lw_path = '/home/tomo/Documents/letswave6'; 
    pica_path = ''; 
    experiment_path = '/DATA2/autocorrelation'; 

else

    error('host %s not set up yet...', hostname);
    
end

fig_path = fullfile(experiment_path, 'figures'); 
data_path = fullfile(experiment_path, 'data'); 
eeg_path = fullfile(experiment_path, 'eeg'); 
resting_eeg_path = fullfile(eeg_path, 'resting_eeg', 'ds004148'); 

%% 

fs = 200; 

pat = [1 0 1 1 1 1 0 1 1 1 0 0]; %  [1 1 1 0 1 1 1 0 1 1 0 0]  [1 0 1 1 1 1 0 1 1 1 0 0]

n_cycles = 20; 

grid_ioi = 0.2; 

trial_dur = n_cycles * length(pat) * grid_ioi; 


%% 

% make sure acf_tools are added because we'll need the get_lag_harmonics
% function from there...
addpath(genpath(acf_tools_path)); 


% min_lag = 0;
% max_lag = (grid_ioi * length(pat)); 
% lags_meter_rel = [0.8, 1.6];
% lags_meter_unrel = [0.6, 1.0, 1.4, 2.0];



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

%% 

max_freq = 5;
max_freq_plot = max_freq + 0.1; 

f0_to_excl = []; 

freq_meter_rel = [1.25 : 1.25 : max_freq]; 
freq_meter_unrel = [1/2.4 : 1/2.4 : max_freq];
freq_meter_unrel = freq_meter_unrel(~ismembertol(freq_meter_unrel, freq_meter_rel, 1e-6)); 

if ~isempty(f0_to_excl)
    freq_meter_rel(mod(freq_meter_rel, f0_to_excl) < ...
        1e4 * eps(min(freq_meter_rel))) = []; 
    freq_meter_unrel(mod(freq_meter_unrel, f0_to_excl) < ...
        1e4 * eps(min(freq_meter_unrel))) = []; 
end

frex = sort([freq_meter_rel, freq_meter_unrel]);

%%

noise_bins = [2, 5]; 

% noies bins for calculating the SNR of the raw spectra
% (harmonic-snippet-zsocre method as used by Rossion)
noise_bins_snr = [3, 13]; 

%%

save_figs = true; 

fontsize = 14; 


%% return structure 

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end


