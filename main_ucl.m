
%%

par = get_par(); 

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath(par.lw_path)); 
addpath(genpath('lib'))

%%

main_snr_vs_nlags('ir_type', 'square'); 
main_snr_vs_nlags('ir_type', 'erp'); 
main_snr_vs_nlags('ir_type', 'erp2'); 

main_only_noise(); 

main_snr_vs_nlags(); 
main_syncrange_eeg_boot(); 
main_syncrange_eeg(); 
main_syncrange_tapping(); 
