function main_fooof_irasa(par, varargin)

parser = inputParser; 

addParameter(parser, 'ir_type', 'square'); % square, erp, erp2
addParameter(parser, 'prepared_noise', []); % square, erp, erp2

parse(parser, varargin{:});

ir_type = parser.Results.ir_type;
noise = parser.Results.prepared_noise;

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%% simulate

noise_exponent = -1.5; 

noise_type = 'eeg'; % eeg, fractal

% number of simulated repetitions 
n_rep = 200; 

if strcmp(noise_type, 'eeg')
    snrs = logspace(log10(0.2), log10(2), 5); 
else
    snrs = logspace(log10(0.2), log10(2), 5); 
end

n_cond = length(snrs); 


%% 

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


%% generate signal

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );

%% genearet noise

if ~isempty(noise)

    n_rep = size(noise, 1); 
    
else
    
    if strcmp(noise_type, 'fractal')

        noise = get_colored_noise2([n_rep, length(x_clean)], par.fs, noise_exponent); 

    elseif strcmp(noise_type, 'eeg')

        noise = prepare_eeg_noise(n_rep, par.trial_dur);    

    else
        error('noise type "%s" not implemented', noise_type);
    end
    
end

%% run

col_names = {'snr', 'only_f0_harmonics', 'r_fooof', 'r_irasa'};

tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);


for only_use_f0_harmonics=[true, false]
    
    rows_all = cell(1, n_cond); 

    parfor i_cond=1:n_cond

        snr = snrs(i_cond); 
        cond_labels{i_cond} = sprintf('%.2g', snr); 

        fprintf('calculating snr %d/%d\n', i_cond, n_cond)

        % scale the noise to the correct SNR 
        x = add_signal_noise(repmat(x_clean, n_rep, 1), noise, snr);

        % get acf
        % -------

        % clean signal
        [acf_clean, lags, ~, mX_clean, freq] = get_acf(x_clean, par.fs);    

        % with aperiodic subtraction using FOOOF
        [acf_fooof, lags, ap_fooof, mX, freq] = get_acf(x, par.fs, ...
                               'rm_ap', true, ...
                               'ap_fit_method', 'fooof', ...
                               'f0_to_ignore', 1/2.4, ...
                               'only_use_f0_harmonics', only_use_f0_harmonics, ...
                               'ap_fit_flims', [0.1, 30]);  

        % with aperiodic subtraction using IRASA
        [acf_irasa, ~, ap_irasa] = get_acf(x, par.fs, ...
                               'rm_ap', true, ...
                               'ap_fit_method', 'irasa', ...
                               'f0_to_ignore', 1/2.4, ...
                               'only_use_f0_harmonics', only_use_f0_harmonics, ...
                               'verbose', 1);  

        % get correelation between ground truth ACF and noise-subtracted one
        r_fooof = correlation(repmat(acf_clean, n_rep, 1), acf_fooof); 
        r_irasa = correlation(repmat(acf_clean, n_rep, 1), acf_irasa); 

        rows_all{i_cond} = [...
            repmat({snr}, n_rep, 1), ...
            repmat({only_use_f0_harmonics}, n_rep, 1), ...
            num2cell(r_fooof), ...
            num2cell(r_irasa) ...
            ];
        
        
        
        figure('color', 'white'); 
        idx = 8; 
        plot(freq, mX(idx, :), 'color', [.5, .5, .5])
        hold on
        plot(freq(2:end), ap_fooof(idx, 2:end), 'b', 'linew', 2)
        plot(freq(2:end), ap_irasa(idx, 2:end), 'r', 'linew', 2)
        box off
        set(gca, 'fontsize', 12); 
        legend({'', 'fooof', 'irasa'}, 'FontSize', 12, 'Box', 'off'); 

    end
    
    % update table 
    rows = vertcat(rows_all{:});

    tbl = [tbl; rows];

end
    
if strcmp(noise_type, 'fractal')
    fname = sprintf('irasaVsFooof_irType-%s_exp-%.1f_nrep-%d', ...
                   ir_type, noise_exponent, n_rep); 
elseif strcmp(noise_type, 'eeg')
    fname = sprintf('irasaVsFooof_irType-%s_noise-eeg_nrep-%d', ...
                   ir_type, n_rep); 
else
    error('noise type "%s" not implemented', noise_type);
end

writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 

% save parameters 
save(fullfile(par.data_path, [fname, '_par.mat']), ...
        'par', 'snrs'); 












