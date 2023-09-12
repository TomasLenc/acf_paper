function match_ylims(all_ax)


lims = [Inf, -Inf]; 
% find
for i_ax=1:length(all_ax)
    ax = all_ax{i_ax}; 
    lims = [min(lims(1), ax.YLim(1)), ...
            max(lims(2), ax.YLim(2))]; 
end

% update
for i_ax=1:length(all_ax)
    ax = all_ax{i_ax}; 
    ax.YLim = lims; 
    ax.YTick = lims; 
end
