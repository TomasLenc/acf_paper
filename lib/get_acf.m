function [acf, lags] = get_acf(data, fs)


% get autocorrelation from FFT
N = size(data, ndims(data)); 
hN = floor( (N + 1) / 2); 
lags = [0 : hN-1] / fs; 

X = fft(data, [], ndims(data)) / N * 2; 
acf = (real(ifft(X .* conj(X), [], ndims(data)))); 

subs_cmd = []; 
subs_cmd.subs = repmat({':'}, 1, ndims(data)); 
subs_cmd.subs{end} = 1:hN;
subs_cmd.type = '()'; 
acf = subsref(acf, subs_cmd); 
