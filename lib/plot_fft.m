function plot_fft(ax, mX, freq, varargin)
     
parser = inputParser; 
addRequired(parser, 'ax'); 
addRequired(parser, 'mX', @isnumeric); 
addRequired(parser, 'freq', @(x) isnumeric(x)); 

addParameter(parser, 'ap', [], @isnumeric); 

addParameter(parser, 'max_freq', max(freq), @isnumeric); 
addParameter(parser, 'freq_meter_rel', [], @isnumeric); 
addParameter(parser, 'freq_meter_unrel', [], @isnumeric); 
addParameter(parser, 'features', [], @isstruct); 
addParameter(parser, 'prec', 100, @isnumeric); 

addParameter(parser, 'col_neutral', repmat(0.5, 1, 3)); 
addParameter(parser, 'col_meter_rel', [0.8706    0.1765    0.1490]); 
addParameter(parser, 'col_meter_unrel', [0.1922    0.5098    0.7412]); 
addParameter(parser, 'col_ap', [0, 0, 0]/255); 
addParameter(parser, 'linew', 2, @isnumeric); 
addParameter(parser, 'linew_ap', 2, @isnumeric); 

parse(parser, ax, mX, freq, varargin{:});

ap = parser.Results.ap; 

max_freq                 = parser.Results.max_freq; 
freq_meter_rel           = parser.Results.freq_meter_rel; 
freq_meter_unrel         = parser.Results.freq_meter_unrel; 
features                 = parser.Results.features; 
prec                     = parser.Results.prec; 

col_neutral     = parser.Results.col_neutral; 
col_meter_rel   = parser.Results.col_meter_rel; 
col_meter_unrel = parser.Results.col_meter_unrel; 
col_ap          = parser.Results.col_ap; 
linew           = parser.Results.linew; 
linew_ap        = parser.Results.linew_ap; 



%%

if isempty(freq_meter_rel) && isempty(freq_meter_unrel)
    y_lims = [min(floor(mX*prec)/prec), max(ceil(mX*prec)/prec)]; 
else
    y_lims = [0, -Inf]; 
end      

hold(ax, 'on');

stem(ax, freq, mX, 'marker', 'none', 'linew', linew, 'color', col_neutral); 

if ~isempty(freq_meter_rel)
    freq_meter_rel_idx = dsearchn(freq', freq_meter_rel'); 
    stem(ax, freq(freq_meter_rel_idx), mX(freq_meter_rel_idx),...
         'marker', 'none', 'linew', linew, 'color', col_meter_rel); 
    y_lims = [0, max(y_lims(2), max(ceil(mX(freq_meter_rel_idx)*prec)/prec))];
end

if ~isempty(freq_meter_unrel)
    freq_meter_unrel_idx = dsearchn(freq', freq_meter_unrel'); 
    stem(ax, freq(freq_meter_unrel_idx), mX(freq_meter_unrel_idx),...
         'marker', 'none', 'linew', linew, 'color', col_meter_unrel); 
    y_lims = [0, max(y_lims(2), max(ceil(mX(freq_meter_unrel_idx)*prec)/prec))];
end

if ~isempty(ap)
    plot(ax, freq, ap, '--', 'marker', 'none', 'linew', linew_ap, 'color', col_ap)
end


x_lims = [0, max_freq]; 
ax.XTick = x_lims; 
ax.XLim = x_lims; 
ax.YTick = y_lims; 
ax.YLim = [floor(y_lims(1)*prec)/prec, ...
           ceil(y_lims(2)*prec)/prec]; 
ax.TickDir = 'out'; 

if ~isempty(features)
    tit = ''; 
    keys = fieldnames(features); 
    for i_key=1:length(keys)
        tit = [tit, sprintf('%s=%.2f   ', keys{i_key}, features.(keys{i_key}))];    
    end
    h_tit = title(ax, tit, 'Interpreter', 'none');
    h_tit.HorizontalAlignment = 'left'; 
    h_tit.Units = 'normalized'; 
    h_tit.Position(1) = 0; 
    h_tit.Position(2) = 1.1; 

end















