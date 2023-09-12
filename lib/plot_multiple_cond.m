function [f, pnl] = plot_multiple_cond(tbl, varargin)

parser = inputParser; 

addParameter(parser, 'ax',  []); 
addParameter(parser, 'colors',  []); 
addParameter(parser, 'plot_legend',  true); 
addParameter(parser, 'ytick_at_means',  false); 
addParameter(parser, 'zero_line',  false); 
addParameter(parser, 'ylim_quantile_cutoff',  0); % between 0 and 1
addParameter(parser, 'prec',  3); % number decimal places 
addParameter(parser, 'point_alpha',  0.4); % opacity of the individual datapoints

addParameter(parser, 'cond_column',  []); 
addParameter(parser, 'feat_column',  []); 
addParameter(parser, 'feat_orig_column',  []); 
addParameter(parser, 'feat_thr_column',  []); 

parse(parser, varargin{:});

ax = parser.Results.ax; 
colors = parser.Results.colors; 
plot_legend = parser.Results.plot_legend; 
ytick_at_means = parser.Results.ytick_at_means; 
zero_line = parser.Results.zero_line; 
prec = parser.Results.prec; 
point_alpha = parser.Results.point_alpha; 

% feat = parser.Results.feat; 
% feat_orig = parser.Results.feat_orig; 
% feat_thr = parser.Results.feat_thr; 
cond_column = parser.Results.cond_column; 
feat_column = parser.Results.feat_column; 
feat_orig_column = parser.Results.feat_orig_column; 
feat_thr_column = parser.Results.feat_thr_column; 

ylim_quantile_cutoff = parser.Results.ylim_quantile_cutoff; 

%% 

conds = unique(tbl{:, cond_column}); 

n_cond = length(conds); 

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

if ~isempty(feat_thr_column)
    
    plot(ax, [0, n_cond+1], [feat_thr.data, feat_thr.data], ...
         '-', 'color', 'red', 'linew', 3)
    ylims(1) = feat_thr.data - abs(feat_thr.data)*0.1;
    ylims(2) = feat_thr.data + abs(feat_thr.data)*0.1;
    
end

if ~isempty(feat_orig_column)
    
    for i_cond=1:n_cond
        
        col = [0.5, 0.5, 0.5];
        
        data = tbl{tbl{:, cond_column} == conds(i_cond), feat_orig_column}; 

        data = unique(data); 
        
        assert(length(data) == 1); 
        
        plot(ax, ...
             [i_cond-0.4, i_cond+0.4], [data, data], ...
             '-', 'color', col, 'linew', 3)
        ylims(1) = min(data - abs(data)*0.1, ylims(1)); 
        ylims(2) = max(data  + abs(data)*0.1, ylims(2)); 
    end
    
end

h = []; 
means = nan(1, n_cond);

if ~isempty(feat_column)
    for i_cond=1:n_cond
        
        data = tbl{tbl{:, cond_column} == conds(i_cond), feat_column}; 
        
        if ~isempty(colors)
            col = colors{i_cond}; 
        else
            col = [120, 20, 166]/255;  
        end
        h(i_cond) = plot_points(ax, i_cond, data, ...
                                'opacity', point_alpha, ...
                                'col', col); 

        ylims(1) = min(quantile(data, ylim_quantile_cutoff), ylims(1)); 
        ylims(2) = max(quantile(data, 1-ylim_quantile_cutoff), ylims(2)); 

        means(i_cond) = mean(data);
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

if plot_legend && ~isempty(feat_column)
    leg = legend(h, cellfun(@num2str, num2cell(conds), 'uni', 0)); 
    leg.Box = 'off'; 
    leg.Position(1) = 0.9278; 
    leg.Position(2) = 0.7700; 
end