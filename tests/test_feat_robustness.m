
% This is quite confusing.... So it seems like when I compute the ACF by
% taking Pearson's correlation coefficient at every single circular lag, it
% gives me the exact same ACF values (exactlty!) as when I take the ACF
% using DFT and then simply dividie the whole ACF by the value at lag 0
% (i.e. variance of x). 

% So that's cool - the ACF we're getting via DFT is just a scaled version
% of the normalized ACF we'd get using Pearson. 

% Now the problem is that ONLY z-score stays invariant to changes in
% unitary response shape. 

% What's crazy is that when I DON'T zscore x before taking the ACF, both
% ratio and contrast work as well - they remain invariant to unitary
% response. 
% 
% However, seems like as soon as I subtract the mean of the
% signal it fucks up the invariance except for zscore, no matter how the
% ACF was calculated and normalized. Even when taken from Pearson's
% correlation, the values don't remain constant across changes in unitary
% response shape.

% So I don't know what's going on in the ACF across the changes in unitary
% response shape... but it cannot be just driven scaling an offset of the
% input signal because we know that Pearson's correlation is invariant to
% that.... I hate not understanding this....!!!

load /Users/tomaslenc/projects_backed_up/acf/data/maxlag-halfTrial_meterRel-0.8_meterUnrel-0.6_1.0_1.4_ignore-0.4/ir-square_ir.mat

par = get_par; 

figure
ax = axes; 

cols = get(ax, 'ColorOrder');

max_lag = 24; 

for i=1:5  
    
    x = data_to_plot(i).x; 
    t = data_to_plot(i).t; 
    fs = 1/data_to_plot(i).t(2); 
                    
    n_lags = max_lag * fs; 
    lags = [0 : n_lags - 1] / fs; 
    acf = nan(1, n_lags); 
    
    for i_lag=1:n_lags
        x1 = x; 
        x2 = circshift(x, -i_lag+1); 
        c = corrcoef(x1, x2); 
        acf(i_lag) = c(1, 2); 

%         acf(i_lag) = dot(x1-mean(x1), x2-mean(x2)) / ...
%                       (std(x1, 1) * std(x2, 1))
        
    end
     
%     figure
    plot(lags, acf); 
    hold on 
    
    [acf, lags] = get_acf(x, fs, ...
                          'normalize_x', true);     
    idx = dsearchn(lags', 2.4); 
    acf = acf ./ acf(idx); 

%     plot(lags, acf); 
%     xlim([0, 10])

    
    feat = get_acf_features(acf, lags, ...
        par.lags_meter_rel,  par.lags_meter_unrel)
end

xlim([0, 10])