function [f, pnl] = plot_multiple_cond(varargin)

parser = inputParser; 

addParameter(parser, 'ax',  []); 
addParameter(parser, 'plot_legend',  true); 
addParameter(parser, 'ytick_at_means',  false); 
addParameter(parser, 'zero_line',  false); 
addParameter(parser, 'ylim_quantile_cutoff',  0); % between 0 and 1
addParameter(parser, 'prec',  3); % number decimal places 
addParameter(parser, 'point_alpha',  0.4); % opacity of the individual datapoints

addParameter(parser, 'feat',  []); 
addParameter(parser, 'feat_orig',  []); 
addParameter(parser, 'feat_thr',  []); 

parse(parser, varargin{:});

ax = parser.Results.ax; 
plot_legend = parser.Results.plot_legend; 
ytick_at_means = parser.Results.ytick_at_means; 
zero_line = parser.Results.zero_line; 
prec = parser.Results.prec; 
point_alpha = parser.Results.point_alpha; 

feat = parser.Results.feat; 
feat_orig = parser.Results.feat_orig; 
feat_thr = parser.Results.feat_thr; 
ylim_quantile_cutoff = parser.Results.ylim_quantile_cutoff; 


n_cond = length(feat); 

ylims = [Inf, -Inf]; 

if isempty(ax)
    f = figure('color', 'white', 'position', [618 374 175 199]); 
    pnl = panel(f); 
    pnl.pack(); 
    pnl.margin = [16, 3, 3, 18];
    pnl.fontsize = 12; 
    ax = pnl(1).select(); 
else
    f = []; 
    pnl = []; 
end

hold(ax, 'on'); 

if ~isempty(feat_thr)
    
    plot(ax, [0, n_cond+1], [feat_thr.data, feat_thr.data], ...
         '-', 'color', 'red', 'linew', 3)
    ylims(1) = feat_thr.data - abs(feat_thr.data)*0.1;
    ylims(2) = feat_thr.data + abs(feat_thr.data)*0.1;
    
end

if length(feat_orig) == 1
    
    plot(ax, [0, n_cond+1], [feat_orig, feat_orig], ...
         '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)
    ylims(1) = feat_orig - abs(feat_orig)*0.1;
    ylims(2) = feat_orig + abs(feat_orig)*0.1;
    
elseif length(feat_orig) == length(feat)
    
    for i_cond=1:n_cond
        % check if we have info about color
        if isfield(feat_orig, 'color')
            col = feat_orig(i_cond).color; 
        else
            col = [0.5, 0.5, 0.5];
        end
        plot(ax, ...
             [i_cond-0.4, i_cond+0.4], ...
             [feat_orig(i_cond).data, feat_orig(i_cond).data], ...
             '-', 'color', col, 'linew', 3)
        ylims(1) = min(feat_orig(i_cond).data - abs(feat_orig(i_cond).data)*0.1, ylims(1)); 
        ylims(2) = max(feat_orig(i_cond).data  + abs(feat_orig(i_cond).data)*0.1, ylims(2)); 
    end
    
end

h = []; 
means = nan(1, n_cond);

if ~isempty(feat)
    for i_cond=1:n_cond
        h(i_cond) = plot_points(ax, i_cond, feat(i_cond).data, ...
                                'opacity', point_alpha, ...
                                'col', feat(i_cond).color); 

        ylims(1) = min(quantile(feat(i_cond).data, ylim_quantile_cutoff), ylims(1)); 
        ylims(2) = max(quantile(feat(i_cond).data, 1-ylim_quantile_cutoff), ylims(2)); 

        means(i_cond) = mean(feat(i_cond).data);
    end

    ax.XLim = [0.5, n_cond+0.5]; 
    if zero_line
        plot(ax.XLim, [0, 0], 'color', [0 0 0], 'linew', 0.5)
    end
    if ylims(1) < ylims(2)
        ax.YLim = ylims; 
    end
    ax.XAxis.Visible = 'off';  
    ax.TickDir = 'out'; 
    if ytick_at_means
        ax.YTick = round(sort(means), prec);
    else
        ax.YTick = [ax.YTick(1), ax.YTick(end)];
    end
    
end

if plot_legend && ~isempty(feat)
    leg = legend(h, {feat.name}); 
    leg.Box = 'off'; 
    leg.Position(1) = 0; 
    leg.Position(2) = 0.8; 
end