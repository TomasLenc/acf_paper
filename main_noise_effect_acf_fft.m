function main_noise_effect_acf_fft(par, varargin)
% clear 
% par = get_par(); 

parser = inputParser; 

addParameter(parser, 'ir_type', 'square'); % square, erp, erp2
addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

ir_type = parser.Results.ir_type;
noise = parser.Results.prepared_noise;


addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))

%% simulate

n_rep = 50; 

emph_levels = 0; % linspace(0, 2, 4)

snrs = logspace(log10(0.01), log10(8), 10); 


%% get IR 

if strcmp(ir_type, 'square')
    ir = get_square_kernel(par.fs, ...
        'duration', par.grid_ioi, ...
        'rampon', 0.010, ...
        'rampoff', 0.010 ...
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


%% find patterns 

n_events = 12; 
n_sounds = 7; 

f_gen_all_pats = @(n, k) dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0'; 
all_possible_pats = f_gen_all_pats(n_events, n_sounds); 
all_possible_pats = flip(all_possible_pats, 1); 

all_good_pats = nan(size(all_possible_pats)); 
z_fft_all_good_pats = nan(size(all_possible_pats, 1), 1); 
z_acf_all_good_pats = nan(size(all_possible_pats, 1), 1); 

c = 1; 
for i_pat=1:size(all_possible_pats, 1)
       
    fprintf('c=%d, pat %d/%d\n', c, i_pat, size(all_possible_pats, 1)); 

    pat = all_possible_pats(i_pat, :); 

    flag = 1; 
    
    for i_shift=0:n_events
        [B, N] = RunLength(circshift(pat, i_shift));
        if any(N(B == 1) > 4) 
            flag = 0; 
        end
    end
    if flag == 0
        continue
    end
    
    for i_shift=0:n_events
        if any(all((all_good_pats - circshift(pat, i_shift)) == 0, 2)) || ...
           any(all((all_good_pats - circshift(flip(pat), i_shift)) == 0, 2))
            flag = 0; 
        end
    end    
    if flag == 0
        continue
    end
    
    [x_clean, t] = get_s(...
                        pat, ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', ir ...
                        );

    [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

    feat_acf_orig = get_acf_features(acf_clean, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

    feat_fft_orig = get_fft_features(mX_clean, freq,...
                            par.freq_meter_rel, par.freq_meter_unrel); 

    z_acf_all_good_pats(c) = feat_acf_orig.z_meter_rel; 
    z_fft_all_good_pats(c) = feat_fft_orig.z_meter_rel; 

    all_good_pats(c, :) = pat; 
    
    c = c+1; 
end

all_good_pats(c:end, :) = []; 
z_acf_all_good_pats(c:end) = []; 
z_fft_all_good_pats(c:end) = []; 


% figure
% [z, idx] = sort(z_acf_all_good_pats); 
% plot(z)
% hold on 
% [z, idx] = sort(z_fft_all_good_pats); 
% plot(z)
% 

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

all_pats = all_good_pats([pat_idx_acf, pat_idx_fft], :); 

%%

if ~isempty(noise)

    n_rep = size(noise, 1); 
    
else
 
    noise = prepare_eeg_noise(n_rep, par.trial_dur);    

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

for i_pat=1:size(all_pats, 1)
    
    % make whole signal 
    [x_clean, t] = get_s(...
                        all_pats(i_pat, :), ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', ir);

    parfor i_snr=1:length(snrs)

        snr = snrs(i_snr); 

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
        [acf_subtracted, ~, ap, ~, ~, par_ap, x_subtr] = get_acf(x, par.fs, ...
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
fname = sprintf('irType-%s_apFitMethod-%s_onlyHarm-%s_nrep-%d_noiseEffectACFvsFFT', ...
                ir_type, ...
                par.ap_fit_method, ...
                jsonencode(par.only_use_f0_harmonics), ...
                n_rep); 
            
writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), 'par', 'all_pats', 'snrs'); 


