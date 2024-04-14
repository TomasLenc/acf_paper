function plot_time(ax, x, t, varargin)

parser = inputParser; 

addParameter(parser, 'col', [120, 20, 166]/255, @isnumeric); 
addParameter(parser, 'xlim', [0, 4.8], @isnumeric); 

parse(parser, varargin{:});

col = parser.Results.col; 
xLims = parser.Results.xlim; 

plot(ax, t, x, 'linew', 2, 'color', col)

ax.XTick = xLims; 
ax.XLim = xLims; 
ax.YTick = ax.YLim; 
ax.TickDir = 'out'; 
