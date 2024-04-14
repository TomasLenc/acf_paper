function h = plot_points(ax, x_coord, vals, varargin)

parser = inputParser; 

addParameter(parser, 'col', [120, 20, 166]/255, @isnumeric); 
addParameter(parser, 'opacity', 0.5, @isnumeric); 

parse(parser, varargin{:});
col = parser.Results.col; 
opacity = parser.Results.opacity; 

hold(ax, 'on')

scatter(ax, repmat(x_coord, 1, length(vals)), vals, ...
        'o', 'MarkerEdgeColor', [0.7, 0.7, 0.7], 'MarkerEdgeAlpha', opacity); 
        
mu = mean(vals);
sd = std(vals); 
ci = sd / sqrt(length(vals)) * norminv(0.95); 
h = plot(ax, x_coord, mu, 'o', 'color', col, 'MarkerFaceColor', col);
plot(ax, [x_coord, x_coord], [mu-sd, mu+sd], 'color', col, 'linew', 2)


