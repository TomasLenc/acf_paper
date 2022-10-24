function plot_time(t, x, cols)

col_time = cols.col_time; 

plot(t, x, 'linew', 2, 'color', col_time)
box off

xLims = [0, 4.8]; 

ax = gca; 
ax.XTick = xLims; 
ax.XLim = xLims; 
ax.YTick = ax.YLim; 

% if strcmp(res.irParams.type,'square')
%     
%     title(sprintf('[%s]  jitter=%.3fs\n duty=%.3fs | rampon=%.3fs | rampoff=%.3fs',...
%                     num2str(res.pattern,'%g '), ...
%                     res.jitter, ...
%                     res.irParams.eventDur, ...
%                     res.irParams.rampon, ...
%                     res.irParams.rampoff)); 
% 
% elseif strcmp(res.irParams.type,'erp')
%     
%     title(sprintf('[%s]  jitter=%.3fs\n  Dur=%.3fs  |  F0 = %s Hz   |   A = %s',...
%                     num2str(res.pattern,'%g '), ...
%                     res.jitter, ...
%                     res.irParams.Dur, ...
%                     num2str(res.irParams.F0), ...
%                     num2str(res.irParams.A))); 
% 
%     
% elseif strcmp(res.irParams.type,'click')
%     
%     title(sprintf('[%s]   jitter=%.3fs\n duration=%.3fs ',...
%                     num2str(res.pattern,'%g '), ...
%                     res.jitter, ...
%                     res.irParams.eventDur)); 
%     
%     
% end