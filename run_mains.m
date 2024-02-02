

par.ir_type = 'square'; 
par.ir = get_ir(par.ir_type, par.fs); 

% main_ir(par);
% 
% main_emph(par);
% 
% main_coch(par);
% 
% main_jitter(par);



%% effect of noise 

par.n_rep = 50; 


% main_snr(par,...
%     'prepared_noise', noise_all_samples);
% 
% main_emph_vs_noise(par, ...
%     'prepared_noise', noise_all_samples);
% 
% 
% 
% main_noiseEffectZscore_ACFvsFFT(par, ...
%     'prepared_noise', noise_all_samples);
% 
% 
% main_noiseEffectDist_ACFvsFFT(par, ...
%     'prepared_noise', noise_all_samples);
% 
% 
% main_noiseEffectDist_allVsOnlyHarm(par, ...
%     'prepared_noise', noise_all_samples); 
% 
main_noiseEffectDist_band(par, ...
    'prepared_noise', noise_all_samples); 
% 
% 
% main_fooof_irasa(par,...
%     'prepared_noise', noise_all_samples); 
% 
% 

% % need loads of samples for this one... 
% par.n_rep = 500; 
% 
% main_only_noise(par,...
%     'prepared_noise', noise_all_samples); 
% 


%% real EEG and tapping data 

% main_lowhigh_eeg(par); 
% 
% main_lowhigh_tap(par); 
% 
% main_infant_eeg(par); 

