function [f, pnl] = plot_multiple_cond(varargin)

parser = inputParser; 

addParameter(parser, 'ax',  []); 
addParameter(parser, 'plot_legend',  true); 
addParameter(parser, 'ylim_quantile_cutoff',  0); % between 0 and 1

addParameter(parser, 'feat',  [], @isstruct); 
addParameter(parser, 'feat_orig',  [], @(x) length(x)==1 | isstruct(x)); 

parse(parser, varargin{:});

ax = parser.Results.ax; 
plot_legend = parser.Results.plot_legend; 

feat = parser.Results.feat; 
feat_orig = parser.Results.feat_orig; 
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

if length(feat_orig) == 1
    
    plot(ax, [0, n_cond+1], [feat_orig, feat_orig], ...
         '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)
    ylims(1) = feat_orig - abs(feat_orig)*0.1;
    ylims(2) = feat_orig + abs(feat_orig)*0.1;
    
elseif length(feat_orig) == length(feat)
    
    for i_cond=1:n_cond
        plot(ax, ...
             [i_cond-0.4, i_cond+0.4], ...
             [feat_orig(i_cond).data, feat_orig(i_cond).data], ...
             '-', 'color', [0.5, 0.5, 0.5], 'linew', 3)
        ylims(1) = min(feat_orig(i_cond).data - abs(feat_orig(i_cond).data)*0.1, ylims(1)); 
        ylims(2) = max(feat_orig(i_cond).data  + abs(feat_orig(i_cond).data)*0.1, ylims(2)); 
    end
    
end

h = []; 
for i_cond=1:n_cond
    h(i_cond) = plot_points(ax, i_cond, feat(i_cond).data, ...
                            'opacity', 0.4, ...
                            'col', feat(i_cond).color); 
                        
    ylims(1) = min(quantile(feat(i_cond).data, ylim_quantile_cutoff), ylims(1)); 
    ylims(2) = max(quantile(feat(i_cond).data, 1-ylim_quantile_cutoff), ylims(2)); 
    
end

ax.XLim = [0.5, n_cond+0.5]; 
ax.YLim = ylims; 
ax.XAxis.Visible = 'off';  
ax.TickDir = 'out'; 
ax.YTick = [ax.YTick(1), ax.YTick(end)];

if plot_legend
    leg = legend(h, {feat.name}); 
    leg.Box = 'off'; 
    leg.Position(1) = 0; 
    leg.Position(2) = 0.8; 
end
