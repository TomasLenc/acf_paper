
if strcmpi(hostname, 'tux')
    
    main03_ir(par, 'ir_type', 'square');
    close all

    main04_snr(par,...
        'prepared_noise', noise_all_samples(1:n_rep_snr, :));
    close all

    main05_emph(par, 'ir_type', 'square');
    close all
    main05_emph(par, 'ir_type', 'erp2');
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
        'prepared_noise', noise_all_samples(1:n_rep_snr_vs_nlags, :)); 

elseif strcmpi(hostname, 'tomo-office-desktop')   
    
    main_only_noise(par,...
        'prepared_noise', noise_all_samples(1:n_rep_only_noise, :)); 
    
    main_syncrange_eeg_boot(par); 
    
    main_syncrange_eeg(par); 
    
    main_syncrange_tapping(par); 
    
    close all
    
end


