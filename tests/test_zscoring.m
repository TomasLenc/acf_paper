
%% 

% Z-score when all beat-unrelated frequecnies are constant. Irrespective of
% the actual value at beat-related frequencies, the z-score at beat
% frequencies will be constant as long as the magnitude at beat frequencies
% is higher than the constant value at beat urelated frequencies (it will
% flip sign if it's lower). 

amps = linspace(1e-3, 1e2, 30); 

x = [5 5 5 5 5 5 5 5 5 5 5 5 5 5]; 

for i=1:length(amps)
    
    x(1) = amps(i); 
    
    z = zscore(x); 
    
    fprintf('amp = %g \t z = %f\n', amps(i), z(1)); 

end


fprintf('\n\n-----------------------------------------------------\n\n'); 


%%

% Problem: the SD of magnitues at beat frequencies will determine the mean
% z-score, even if the mean magnitude at beat frequencies remains constant.
% This means that the measure is not invariant to the shape of the periodic
% signal - some periodic signals will have wider/narrower distribution of
% energy across harmonics (e.g. isochronous dirac implses have 0 SD, but
% erp-like response has much larger SD). 

% However, this should not be a problem in ACF - as long as the unitary
% response doesn't last longer than grid interval, the ACF values at
% beat-related/unrelated lags will remain invariant to the shape of the
% unitary response. 

idx_meter_rel = [1 2 3]; 

x = [6 6 6 5 5 5]; 
z = zscore(x); 
mean(x(idx_meter_rel))
mean(z(idx_meter_rel))

x = [5 6 7 5 5 5]; 
z = zscore(x); 
mean(x(idx_meter_rel))
mean(z(idx_meter_rel))


%%

% Still, one thing I'm trying to understand is whether or not we're
% discriminating against some rhythms that may have comparable periodic
% recurrence to other rhythms, but based on their structure, the SD of acf
% across beat-related and beat-unrelated lags is just different. Can't
% think of good examples and I don't have clear intuition in my head right
% now. 




%% 

% Problem: when all values are near constant, noise will have a huge
% effect. 




% It is the same with ratio, percent difference, or Michelson contrast 




