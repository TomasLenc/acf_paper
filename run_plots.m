
clear

%%

% 
% sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4'; 
% 
% par = get_par(); 
% par.data_path = fullfile(par.data_path, sel_name); 
% 
% par.ir_type = 'square'; 
% plot_ir(par); 
% par.ir_type = 'erp'; 
% plot_ir(par); 
% 
% par.ir_type = 'square'; 
% plot_emph(par); 
% par.ir_type = 'erp2'; 
% plot_emph(par); 
% 
% par.ir_type = 'square'; 
% plot_snr(par); 
% par.ir_type = 'erp2'; 
% plot_snr(par); 
% 
% par.ir_type = 'square'; 
% plot_jitter(par); 
% par.ir_type = 'erp2'; 
% plot_jitter(par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-frontocentral_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_lowhigh(par_eeg.par); 
% plot_lowhigh(par_eeg.par, 'chunk');
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_infant(par_eeg.par); 
% 
% %%
% 
% sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 
% 
% par = get_par(); 
% par.data_path = fullfile(par.data_path, sel_name); 
% 
% par.ir_type = 'square'; 
% plot_ir(par); 
% par.ir_type = 'erp'; 
% plot_ir(par); 
% 
% par.ir_type = 'square'; 
% plot_emph(par); 
% par.ir_type = 'erp2'; 
% plot_emph(par); 
% 
% par.ir_type = 'square'; 
% plot_snr(par); 
% par.ir_type = 'erp2'; 
% plot_snr(par); 
% 
% par.ir_type = 'square'; 
% plot_jitter(par); 
% par.ir_type = 'erp2'; 
% plot_jitter(par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-frontocentral_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_lowhigh(par_eeg.par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_infant(par_eeg.par); 
% 
% 
% 
% %%
% 
% sel_name = 'maxlag-halfTrial_meterRel-0.4_meterUnrel-0.6_1.0_1.4'; 
% 
% par = get_par(); 
% par.data_path = fullfile(par.data_path, sel_name); 
% 
% par.ir_type = 'square'; 
% plot_ir(par); 
% par.ir_type = 'erp'; 
% plot_ir(par); 
% 
% par.ir_type = 'square'; 
% plot_emph(par); 
% par.ir_type = 'erp2'; 
% plot_emph(par); 
% 
% par.ir_type = 'square'; 
% plot_snr(par); 
% par.ir_type = 'erp2'; 
% plot_snr(par); 
% 
% par.ir_type = 'square'; 
% plot_jitter(par); 
% par.ir_type = 'erp2'; 
% plot_jitter(par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-frontocentral_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_lowhigh(par_eeg.par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_infant(par_eeg.par); 
% 
% 
% %%
% 
% 
% sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.2'; 
% 
% par = get_par(); 
% par.data_path = fullfile(par.data_path, sel_name); 
% 
% par.ir_type = 'square'; 
% plot_ir(par); 
% par.ir_type = 'erp'; 
% plot_ir(par); 
% 
% par.ir_type = 'square'; 
% plot_emph(par); 
% par.ir_type = 'erp2'; 
% plot_emph(par); 
% 
% par.ir_type = 'square'; 
% plot_snr(par); 
% par.ir_type = 'erp2'; 
% plot_snr(par); 
% 
% par.ir_type = 'square'; 
% plot_jitter(par); 
% par.ir_type = 'erp2'; 
% plot_jitter(par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-frontocentral_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_lowhigh(par_eeg.par); 
% 
% 
% par_eeg = load(fullfile(par.data_path, ...
%     'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
% par_eeg.par.data_path = fullfile(par.data_path); 
% 
% plot_infant(par_eeg.par); 
% 



%%


sel_name = 'maxlag-11lags_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

par.ir_type = 'square'; 
plot_ir(par); 

par.ir_type = 'square'; 
plot_emph(par); 

par.ir_type = 'square'; 
plot_snr(par); 

par.ir_type = 'square'; 
plot_jitter(par); 


par_eeg = load(fullfile(par.data_path, ...
    'exp-lowhigh_apFitMethod-irasa_onlyHarm-true_roi-frontocentral_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_lowhigh(par_eeg.par); 


par_eeg = load(fullfile(par.data_path, ...
    'exp-infant_apFitMethod-irasa_onlyHarm-true_roi-front_eegIndividual_par')); 
par_eeg.par.data_path = fullfile(par.data_path); 

plot_infant(par_eeg.par); 



