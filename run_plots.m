
clear

%%


sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

plot_ir(par); 

plot_emph(par); 

plot_snr(par); 

plot_jitter(par); 


%%

sel_name = 'maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4'; 

par = get_par(); 
par.data_path = fullfile(par.data_path, sel_name); 

plot_ir(par); 

plot_emph(par); 

plot_snr(par); 

plot_jitter(par); 
