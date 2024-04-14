
clear



%%

only_use_f0_harmonics   = [false, true, true]; 
keep_band               = [false, false, true]; 

sel_name_templ = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-%s_keepBand-%s'; 

for i_sel=1:length(only_use_f0_harmonics)

    sel_name = sprintf(sel_name_templ, ...
                       jsonencode(only_use_f0_harmonics(i_sel)), ...
                       jsonencode(keep_band(i_sel))); 
    
    par = get_par(); 
    
    par.data_path = fullfile(par.data_path, sel_name); 

    par.only_use_f0_harmonics = only_use_f0_harmonics(i_sel); 
    
    if keep_band(i_sel)
        par.ap_band_around_harmonics = [2, 5]; 
    else
        par.ap_band_around_harmonics = [1, 1];         
    end
    
    % ----------------------------
    
    par.ir_type = 'square'; 
    plot_ir(par); 
    par.ir_type = 'erp'; 
    plot_ir(par); 

    par.ir_type = 'square'; 
    plot_emph(par); 
    par.ir_type = 'erp2'; 
    plot_emph(par); 

    par.ir_type = 'square'; 
    plot_jitter(par); 

    % lowhigh roi-front
    par_eeg = load(fullfile(par.data_path, ...
        sprintf('exp-lowhigh_apFitMethod-irasa_onlyHarm-%s_roi-front_eegIndividual_par', ...
                jsonencode(only_use_f0_harmonics(i_sel))))); 
    par_eeg.par.data_path = fullfile(par.data_path); 

    plot_lowhigh(par_eeg.par); 

    % infant roi-front
    par_eeg = load(fullfile(par.data_path, ...
        sprintf('exp-infant_apFitMethod-irasa_onlyHarm-%s_roi-front_eegIndividual_par', ...
                jsonencode(only_use_f0_harmonics(i_sel))))); 
    par_eeg.par.data_path = fullfile(par.data_path); 

    plot_infant(par_eeg.par); 

    % plot example of nonrepeating sequence 
    res = load(fullfile(par.data_path, 'ir-square_ir.mat')); 
    res.par.max_lag = 12; 
    f = plot_nonrep_seq(res.par); 
    save_fig(f, fullfile(par.data_path, 'ir-square_nonrep.svg'));
    
    close all

end








