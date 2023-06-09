function main_syncrange_eeg(par)
% clear 
% par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))


%% parameters

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.05; 

%%

tbl_subj = readtable('xpsyncrange_preproc_dirnames.txt'); 

rhythms = {'31', '26', '19', '37', '42', '6', '17', '41'}; 

cond_type = 'rhythm'; 

n_rhythms = length(rhythms); 

% colors
cmap_name = 'Set1'; 
colors = num2cell(brewermap(n_rhythms, cmap_name), 2); 

%% allocate table

col_names = {
    'subject', 'rhythm', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_sound', 'z_meter_acf_sound', ...
    'z_snr', 'ap_offset', 'ap_exponent' ...
    };

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


%% RUN


%% allocate

feat_acf_s = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
    'contrast_meter_rel', []); 

feat_acf = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
    'contrast_meter_rel', []); 

feat_acf_subtracted = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
    'contrast_meter_rel', []); 

feat_fft_s = struct('z_meter_rel', []); 

feat_fft = struct('z_meter_rel', []); 

feat_fft_subtracted = struct('z_meter_rel', []); 

feat_ap = struct('offset', [], 'exponent', []); 

%% run

for i_rhythm=1:n_rhythms
    
    rhythm_id = rhythms{i_rhythm};
    
    fprintf('processing rhythm %s\n', rhythm_id);


    %% load data
    
    fpath_eeg = fullfile(par.eeg_path, 'syncrange'); 
    
    eeg = load(fullfile(fpath_eeg, ...
        sprintf('exp-syncrange_rhythm-%s_nTrials-9_eeg.mat', ...
                rhythm_id)));
    
    data = eeg.data;
    fs = eeg.fs;
    t = [0 : size(eeg.data, 2) - 1] / fs; 
        
    %% load stimulus
    
    fpath_stim = '/DATA1/XPSyncRange/ptb/eeg_7set/to_load'; 
    
    d = dir(fullfile(fpath_stim, sprintf('*rhythm%s_*.wav', rhythm_id)));
    
    [s, fs_s] = audioread(fullfile(fpath_stim, d.name));
    s = s(:, 1)'; 
    
    env = abs(hilbert(s)); 
    
    t_s = [0 : length(s) - 1] / fs_s; 
        
    %% process stimulus
    
    % get ACF (withuout aperiodic subtraction)
    [acf_s, lags_s, ~, mX_s, freq_s] = get_acf(env, fs_s);    
                           
    % get ACF features
    feat_acf_s(i_rhythm) = get_acf_features(...
                                acf_s, lags_s, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    
    
    % get features for the raw spectra                                    
    feat_fft_s(i_rhythm) = get_fft_features(mX_s, freq_s, ...
                                     par.freq_meter_rel, par.freq_meter_unrel); 
    
    
    %% process EEG
    
    n_sub = size(data, 1); 
            
    % get acf
    % -------
                                       
    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(data, fs);    
                                   
    mX_subtracted = subtract_noise_bins(mX, par.noise_bins(1),  par.noise_bins(2)); 
    
    % with aperiodic subtraction    
    [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr, optim_exitflag] = ...
                                get_acf(data, fs, ...
                                       'rm_ap', true, ...
                                       'f0_to_ignore', 1 / 2.4, ...
                                       'min_freq', 0.1, ...
                                       'max_freq', 9);      
    if any(~optim_exitflag)
        warning('ap-fit didnt converge %d/%d reps', sum(~optim_exitflag), n_rhythms); 
    end
    
    feat_ap(i_rhythm).offset = cellfun(@(x) x(1), par_ap);            
    feat_ap(i_rhythm).exponent = cellfun(@(x) x(2), par_ap);            
                                   
    % get features
    % ------------

    feat_acf(i_rhythm) = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    

    feat_acf_subtracted(i_rhythm) = get_acf_features(acf_subtracted, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel); 

    % get features for the raw spectra                                    
    tmp = get_fft_features(mX, freq, par.freq_meter_rel, par.freq_meter_unrel); 
    feat_fft(i_rhythm).z_meter_rel = tmp.z_meter_rel; 
                                        
    feat_fft(i_rhythm).z_snr = get_z_snr(mX, freq, par.frex, ...
                                       par.noise_bins_snr(1), ...
                                       par.noise_bins_snr(2)); 

    % get features for the 1/f-subtracted spectra                                    
    feat_fft_subtracted(i_rhythm) = get_fft_features(mX_subtracted, freq, ...
                                           par.freq_meter_rel, par.freq_meter_unrel);
    
                             
    % plot example 
    % ------------
    
    if i_rhythm==1
        f = figure('color','white', ...
                   'position', [95, 67, 1062, 170 * n_rhythms]); 
        pnl_example = panel(f); 
        pnl_example.pack('v', n_rhythms); 
        pnl_example.margin = [5, 10, 25, 25]; 
    end
    rep_to_plot_idx = 9; 
    plot_example(data(rep_to_plot_idx, :), t, ...
                     acf(rep_to_plot_idx, :), lags, ...
                     ap(rep_to_plot_idx, :), ...
                     mX(rep_to_plot_idx, :), freq, ...
                     par.lags_meter_rel, par.lags_meter_unrel, ...
                     par.freq_meter_rel, par.freq_meter_unrel, ...
                     'pnl', pnl_example(i_rhythm), ...
                     'subplot_proportions', [50, 17, 33], ...
                     'max_t', 40.8, ...
                     'min_lag', 0.2, ...
                     'max_lag', par.max_lag, ...
                     'max_freq', par.max_freq_plot, ...
                     'plot_time_xaxis', i_rhythm == n_rhythms, ...
                     'plot_xlabels', i_rhythm == n_rhythms, ...
                     'plot_xticks', i_rhythm == n_rhythms, ...
                     'plot_features', false, ...
                     'mX_subtr', mX_subtracted(rep_to_plot_idx, :), ...
                     'acf_subtr', acf_subtracted(rep_to_plot_idx, :), ...
                     'time_col', colors{i_rhythm}, ...
                     'prec', 1e9, ...
                     'fontsize', par.fontsize, ...
                     'normalize_acf_for_plotting', false);                                        
    f.Name = rhythms{i_rhythm};     
    pnl_example(i_rhythm).margintop = 25; 

end
    
fname = sprintf('exp-syncrange_response-eegIndividual_examples.svg'); 
save_fig(f, fullfile(par.data_path, fname))

%% save table

fname = sprintf('exp-syncrange_eegIndividual'); 
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 




