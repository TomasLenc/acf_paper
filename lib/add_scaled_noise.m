function eeg = add_scaled_noise(eeg, SNR, fs, noise_exp)
% INPUT: 
%     eeg             matrix, size m (trials) x n (timepoints) 


% get rms of EEG 
eeg_rms = rms(eeg,2); 
N = length(eeg); 

for triali=1:size(eeg,1)
    % generate pink noise
    noise = get_colored_noise(N, fs, noise_exp); 
    % scale it to achieve requested SNR
    noise = (noise / rms(noise)) * (eeg_rms(triali)/SNR); 
    % add noise to EEG signal
    eeg(triali,:) = eeg(triali,:) + noise;     
end

