function plot_acf(acf, ap, lags_time, lags_meter_rel, lags_meter_unrel, ...
                  z_meter_rel, ratio_meter_rel, contrast_meter_rel, ...
                  prec, cols)
     
% % get autocorrelation from matlab builtin function
% [acf, lags] = xcorr(res.s); 
% lagsTime = lags/res.fs; 

lags_meter_rel_idx = dsearchn(lags_time', lags_meter_rel'); 
lags_meter_unrel_idx = dsearchn(lags_time', lags_meter_unrel'); 
 
y_lims = [min(floor(acf*prec)/prec), ...
          max(ceil(acf*prec)/prec)]; 

plot(lags_time, acf, 'marker', 'none', 'linew', 2, 'color', cols.col_acf)
hold on
if ~isempty(ap)
    plot(lags_time, ap, '--', 'marker', 'none', 'linew', 2, 'color', cols.col_ap)
end

hold on
plot([lags_time(lags_meter_rel_idx); lags_time(lags_meter_rel_idx)], y_lims,...
        'r--', 'linew', 2, 'color', cols.col_meter_rel)
plot([lags_time(lags_meter_unrel_idx); lags_time(lags_meter_unrel_idx)], y_lims,...
        'r--', 'linew', 2, 'color', cols.col_meter_unrel)
    
ax = gca; 
x_lims = [0, lags_time(end)]; 
ax.XTick = x_lims; 
ax.XLim = x_lims; 
ax.YTick = y_lims; 
ax.YLim = [floor(y_lims(1)*prec)/prec, ...
           ceil(y_lims(2)*prec)/prec]; 

tit = ''; 
if ~isempty(z_meter_rel)
    tit = [tit, sprintf('z = %.3f    ', z_meter_rel)]; 
end
if ~isempty(ratio_meter_rel)
    tit = [tit, sprintf('ratio = %.3f    ', ratio_meter_rel)]; 
end
if ~isempty(contrast_meter_rel)
    tit = [tit, sprintf('contr = %.3f    ', contrast_meter_rel)]; 
end
title(tit)















