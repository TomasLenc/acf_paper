
% run many simulations where the shape of the unitary event is generated
% randomly, and the features from acf are extracted to see which measure is
% invariant to uniary response shape
main_ir_rand(par);

% -----------------------------------------------
% analyses using 'erp' shaped kernel 

par.ir_type = 'erp'; 
par.ir = get_ir(par.ir_type, par.fs); 

main_ir(par);

% -----------------------------------------------
% analyses using 'erp2' shaped kernel 

par.ir_type = 'erp2'; 
par.ir = get_ir(par.ir_type, par.fs); 

main_emph(par); 
main_jitter(par); 

% -----------------------------------------------
% analyses using 'square wave' shaped kernel 

par.ir_type = 'square'; 
par.ir = get_ir(par.ir_type, par.fs); 

main_ir(par);
main_emph(par); 
main_jitter(par); 

% -----------------------------------------------
% effect of noise

par.n_rep = 50; 

par.ir_type = 'square'; 
par.ir = get_ir(par.ir_type, par.fs); 

% analyse effect of noise on z-score at beat-rel freqs/lags
main_noiseEffectZscore_ACFvsFFT(par, ...
    'prepared_noise', noise_all_samples);

% analyse effect of noise on distance from ground truth (prepare data for
% comparison of FFT and ACF sensitivity to noise level)
main_noiseEffectDist_ACFvsFFT(par, ...
    'prepared_noise', noise_all_samples);

% analyse effect of noise on distance from ground truth (compare the
% efect of noise with/without zeroing-out step in the ACF cleaning
% procedure)
main_noiseEffectDist_allVsOnlyHarm(par, ...
    'prepared_noise', noise_all_samples); 

% analyse effect of noise on distance from ground truth (effect of leaving
% a narrow band around each response frequency during the zero-out noise
% correction step)
main_noiseEffectDist_band(par, ...
    'prepared_noise', noise_all_samples); 

% analyse effect of noise on distance from ground truth (compare the
% performance of FOOOF vs. IRASA in correcting for 1/f)
main_noiseEffectDist_FooofVsIrasa(par,...
    'prepared_noise', noise_all_samples); 


%% pure noise

% need loads of samples for this one... 
par.n_rep = 500; 

par.ir_type = 'square'; 
par.ir = get_ir(par.ir_type, par.fs); 

main_only_noise(par,...
    'prepared_noise', noise_all_samples); 


%% cochlear model 

main_coch(par);


%% real EEG and tapping data 

main_lowhigh_eeg(par); 

main_lowhigh_tap(par); 

main_infant_eeg(par); 

