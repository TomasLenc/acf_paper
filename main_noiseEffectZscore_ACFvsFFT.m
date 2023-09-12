function main_noiseEffectZscore_ACFvsFFT(par, varargin)
% This code first finds generates patterns separately for FFT and ACF, such
% that their resulting meter zscore (under the selected impulse response and
% frex/lags of interest) is near pre-defined values. Then it simulates how the
% meter zscore is affected by noise. The goal is to show how adding more and
% more noise pulls the estimate to zero. 

parser = inputParser; 

addParameter(parser, 'prepared_noise', []); 

parse(parser, varargin{:});

noise = parser.Results.prepared_noise;


%% simulate

par.snrs = logspace(log10(0.01), log10(8), 10); 


%% find patterns 

n_events = 12; 
n_sounds = 7; 
max_group_size = 4; 

all_good_pats = find_all_patterns(n_events, n_sounds, max_group_size); 

z_fft_all_good_pats = nan(size(all_good_pats, 1), 1); 
z_acf_all_good_pats = nan(size(all_good_pats, 1), 1); 

for i_pat=1:size(all_good_pats, 1)
       
    pat = all_good_pats(i_pat, :); 
    
    [x_clean, t] = get_s(...
                        pat, ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', par.ir ...
                        );

    [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

    feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

    feat_fft_orig = get_fft_features(mX_clean, freq,...
                            par.freq_meter_rel, par.freq_meter_unrel); 

    z_acf_all_good_pats(i_pat) = feat_acf_orig.z_meter_rel; 
    z_fft_all_good_pats(i_pat) = feat_fft_orig.z_meter_rel; 

end

% figure
% [z, idx] = sort(z_acf_all_good_pats); 
% plot(z)
% hold on 
% [z, idx] = sort(z_fft_all_good_pats); 
% plot(z)
% 

% idx = find(abs(z_acf_all_good_pats - z_fft_all_good_pats) < 0.3); 
% maxlim = 2.5; 
% minlim = -maxlim; 
% figure
% plot(z_acf_all_good_pats, z_fft_all_good_pats, 'ko'); 
% hold on
% plot(z_acf_all_good_pats(idx), z_fft_all_good_pats(idx),'ro'); 
% plot([minlim, maxlim], [0 0], 'k', 'linew', 1); 
% plot([0 0], [minlim, maxlim], 'k', 'linew', 1); 
% plot([minlim, maxlim], [minlim, maxlim], 'k:', 'linew', 2); 
% axis square
% xlabel('z meter acf'); 
% ylabel('z meter fft'); 
% xticks([minlim : 0.5 : maxlim]); 
% yticks([minlim : 0.5 : maxlim]); 
% xlim([minlim, maxlim])
% ylim([minlim, maxlim])


target_z = [-0.5, 0, 0.5]; 

n_pat_per_z = 3; 

pat_idx_acf = []; 
pat_idx_fft = []; 
for i=1:length(target_z)
    [err, idx] = sort(abs(z_acf_all_good_pats - target_z(i))); 
    pat_idx_acf = [pat_idx_acf; idx(1:n_pat_per_z)]; 
    
    [err, idx] = sort(abs(z_fft_all_good_pats - target_z(i))); 
    pat_idx_fft = [pat_idx_fft; idx(1:n_pat_per_z)]; 
end

z_acf_all_good_pats(pat_idx_acf)
z_fft_all_good_pats(pat_idx_fft)

%%

n_rep = par.n_rep; 

par.all_pats = all_good_pats([pat_idx_acf, pat_idx_fft], :); 


%% prepare noise
if size(noise, 1) < n_rep
    error('you requested %d samples but provided only noise for %s...', ...
          n_rep, size(noise, 1)); 
else
    noise = noise(1:n_rep, :); 
end

%%

col_names = {
    'pat', 'selected_for', 'snr', 'sample', ...
    'z_meter_fft_raw', 'z_meter_acf_raw', ...
    'z_meter_fft_subtr', 'z_meter_acf_subtr', ...
    'z_meter_fft_orig', 'z_meter_acf_orig', ...
    'z_snr' ...
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

        snr = par.snrs(i_snr); 

        fprintf('pat %d, snr %.1f\n', i_pat, snr); 

        % scale the noise to the correct SNR 
        x = add_signal_noise(x_clean, noise, snr);

        % get acf
        % -------

        % clean signal
        [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

        % withuout aperiodic subtraction    
        [acf, ~, ~, mX, ~] = get_acf(x, par.fs);    

        mX_subtracted = subtract_noise_bins(mX, ...
                                par.noise_bins(1), par.noise_bins(2)); 

        % with aperiodic subtraction    
        [acf_subtracted] = get_acf(x, par.fs, ...
                       'rm_ap', true, ...
                       'ap_fit_method', par.ap_fit_method, ...
                       'f0_to_ignore', par.f0_to_ignore, ...
                       'only_use_f0_harmonics', par.only_use_f0_harmonics, ...
                       'ap_fit_flims', par.ap_fit_flims ...
                       );

        % get features
        % ------------

        feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel);         

        feat_acf = get_acf_features(acf, lags, ...
                                    par.lags_meter_rel, par.lags_meter_unrel);    

        feat_acf_subtracted = get_acf_features(acf_subtracted, lags, ...
                                     par.lags_meter_rel, par.lags_meter_unrel); 


        % get features for the clean spectra
        feat_fft_orig = get_fft_features(mX_clean, freq,...
                                par.freq_meter_rel, par.freq_meter_unrel); 

        % get features for the raw spectra                                    
        tmp = get_fft_features(mX, freq, ...
                        par.freq_meter_rel, par.freq_meter_unrel); 
        feat_fft = []; 
        feat_fft.z_meter_rel = tmp.z_meter_rel; 

        feat_fft.z_snr = get_z_snr(mX, freq, par.frex, ...
                                           par.noise_bins_snr(1), ...
                                           par.noise_bins_snr(2)); 

        % get features for the 1/f-subtracted spectra                                    
        feat_fft_subtracted = get_fft_features(mX_subtracted, freq, ...
                                par.freq_meter_rel, par.freq_meter_unrel);

        if i_pat <= n_pat_per_z * length(target_z)
            selected_for = 'acf'; 
        else
            selected_for = 'fft'; 
        end
        rows = [...
            repmat({i_pat}, n_rep, 1), ...
            repmat({selected_for}, n_rep, 1), ...
            repmat({snr}, n_rep, 1), ...
            num2cell([1:n_rep]'), ...
            num2cell(feat_fft.z_meter_rel), ...
            num2cell(feat_acf.z_meter_rel), ...
            num2cell(feat_fft_subtracted.z_meter_rel), ...
            num2cell(feat_acf_subtracted.z_meter_rel), ...
            repmat({feat_fft_orig.z_meter_rel}, n_rep, 1), ...
            repmat({feat_acf_orig.z_meter_rel}, n_rep, 1), ...
            num2cell(feat_fft.z_snr)...
            ];

        tbl = [tbl; rows]; 

    end
end

%%

% save table
fname = sprintf('irType-%s_apFitMethod-%s_onlyHarm-%s_noiseEffectZscoreACFvsFFT', ...
                par.ir_type, ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics)); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par'); 


