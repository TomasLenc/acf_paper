
clear

%%

sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-false_keepBand-false'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = false; 
par.ap_band_around_harmonics = [1, 1]; 

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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-false_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-false_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all


%%

sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-false'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [1, 1]; 

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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all

%%

sel_name = 'maxlag-2.4_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-true'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [2, 5]; 


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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all

%%
%%
%%




%%

sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-false_keepBand-false'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = false; 
par.ap_band_around_harmonics = [1, 1]; 

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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-false_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-false_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all


% plot example of nonrepeating sequence 
res = load(fullfile(par.data_path, 'ir-square_ir.mat')); 
res.par.max_lag = 12; 
f = explain_nonrep_seq(res.par); 
save_fig(f, fullfile(par.data_path, 'ir-square_nonrep.svg'));


%%

sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-false'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [1, 1]; 

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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all

%%

sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4_zeroOut-true_keepBand-true'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.only_use_f0_harmonics = true; 
par.ap_band_around_harmonics = [2, 5]; 


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
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 

% infant roi-front
par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 
 
close all








