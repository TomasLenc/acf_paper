function par = get_par()

acf_tools_path = '/datadisk/projects_git_dl/acf_tools/src'; 
rnb_tools_path = '/datadisk/projects_git_dl/rnb_tools/src'; 

experiment_path = '/datadisk/projects_backed_up/autocorrelation'; 
fig_path = fullfile(experiment_path, 'figures'); 
data_path = fullfile(experiment_path, 'data'); 

lw_path = '/datadisk/projects_backed_up/autocorrelation/lib_external/letswave6'; 
pica_path = '/datadisk/projects_backed_up/autocorrelation/lib_external/piCA'; 

%% 

fs = 200; 

pat = [1 0 1 1 1 1 0 1 1 1 0 0]; %  [1 1 1 0 1 1 1 0 1 1 0 0]  [1 0 1 1 1 1 0 1 1 1 0 0]

n_cycles = 16; 

grid_ioi = 0.2; 

%% 

% autocorrelation lags (in seconds) that are considered meter-related and
% meter-unrelated
min_lag = 0;
max_lag = (grid_ioi * length(pat) * n_cycles) / 2; 

% meter-related lags 
% ------------------

% Make sure there's no overlap with muiltiples of meter-unrelated lags, and 
% also the pattern repetition period. 
lags_meter_rel = get_lag_harmonics(...
                            0.8, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.6, 1.0, 2.4] ...
                            ); 

% meter-unrelated lags 
% --------------------

% Make sure there's no overlap with muiltiples of meter-related lags 
% (even 0.4 seconds!), and also the pattern repetition period. 

lags_meter_unrel_1 = get_lag_harmonics(...
                            0.6, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.4, 2.4] ...
                            ); 

lags_meter_unrel_2 = get_lag_harmonics(...
                            1.0, ...
                            max_lag,...
                            'lag_harm_to_exclude', [0.4, 2.4] ...
                            ); 

lags_meter_unrel = uniquetol([lags_meter_unrel_1, lags_meter_unrel_2], 1e-8); 


% make sure one more time that there's no overlap between meter-rel and -unrel !!! 
assert(~any( min(abs(bsxfun(@minus, lags_meter_rel', lags_meter_unrel))) < 1e-9 ))


% you can separately set meter-unrelated lags on the left and right (this is
% used when checking for spurious results)
lags_meter_unrel_left = [0.6]; 
lags_meter_unrel_right = [1.0]; 

%% 

freq_meter_rel = [1.25 : 1.25 : 5]; 
freq_meter_unrel = setdiff(1/2.4 * [1:12], freq_meter_rel); 

frex = [freq_meter_rel, freq_meter_unrel];

max_freq_plot = 5.5; 
noise_bins = [2, 5]; 

% noies bins for calculating the SNR of the raw spectra
% (harmonic-snippet-zsocre method as used by Rossion)
noise_bins_snr = [3, 13]; 

%%

fontsize = 14; 


%% return structure 

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end


