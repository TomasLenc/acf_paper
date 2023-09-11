
if strcmpi(hostname, 'tux')
    
    main03_ir(par, 'ir_type', 'square');
    close all

    main04_snr(par,...
        'ir_type', 'square', ...
        'prepared_noise', noise_all_samples(1:n_rep_snr, :));
    close all
    main04_snr(par,...
        'ir_type', 'erp2', ...
        'prepared_noise', noise_all_samples(1:n_rep_snr, :));
    close all

    main05_emph(par, 'ir_type', 'square');
    close all
    main05_emph(par, 'ir_type', 'erp2');
    close all
    
    main08_coch(par);
    close all
    
    main10_jitter(par);
    close all

    main_emph_vs_noise(par, ...
        'ir_type', 'square', ...
        'prepared_noise', noise_all_samples(1:n_rep_emph_vs_noise, :));
    
    main_emph_vs_noise(par, ...
        'ir_type', 'erp2', ...
        'prepared_noise', noise_all_samples(1:n_rep_emph_vs_noise, :));

    main_snr_vs_nlags(par,...
        'ir_type', 'square', ...
        'prepared_noise', noise_all_samples(1:n_rep_snr_vs_nlags, :)); 

    main_snr_vs_nlags(par,...
        'ir_type', 'erp2', ...
        'prepared_noise', noise_all_samples(1:n_rep_snr_vs_nlags, :)); 

    main_fooof_irasa(par,...
        'ir_type', 'square', ...
        'prepared_noise', noise_all_samples(1:100, :)); 
    
    
elseif strcmpi(hostname, 'tomo-office-desktop')   
    
%     main_only_noise(par,...
%         'prepared_noise', noise_all_samples(1:n_rep_only_noise, :)); 
    
    par.ap_fit_method = 'irasa'; 
    main_noise_effect_acf_fft(par, ...
        'ir_type', 'square', ...
        'prepared_noise', noise_all_samples(1:50, :));
    

%     % update maximum lag of interest according to eeg trial duration (40.8 s)
%     par.max_lag = 20.4; 
%     par.lags_meter_rel = par.lags_meter_rel(par.lags_meter_rel < par.max_lag); 
%     par.lags_meter_unrel = par.lags_meter_unrel(par.lags_meter_unrel < par.max_lag); 
%     
%     main_syncrange_eeg_boot(par); 
%     
%     main_syncrange_eeg(par); 
%     
%     main_syncrange_tapping(par); 
    
    close all
    
end


