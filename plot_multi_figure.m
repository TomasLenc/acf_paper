function [f, pnl] = plot_multi_figure(data_path, fname, cmap_name, cond_colname, feat_to_plot, varargin)

parser = inputParser; 

addParameter(parser, 'plot_subtr',  true); 
addParameter(parser, 'varargin_for_examples',  {}); 
addParameter(parser, 'varargin_for_points',  {}); 

parse(parser, varargin{:});

plot_subtr = parser.Results.plot_subtr; 
varargin_for_examples = parser.Results.varargin_for_examples; 
varargin_for_points = parser.Results.varargin_for_points; 

%% load 

res = load(fullfile(data_path, [fname, '.mat'])); 
data_to_plot = res.data_to_plot; 
par = res.par; 

tbl = readtable(fullfile(data_path, [fname, '.csv'])); 

n_cond = length(data_to_plot); 

%%

% colors
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(end-n_cond+1:end, :); 

%%

if par.max_lag <= 5
    fig_size = [146 1340 900 700]; 
    example_subplot_proportions = [60, 15, 25];
    min_n_pack_feat_pnl = 4; 
else
    fig_size = [146 1340 1329 700]; 
    example_subplot_proportions = [35, 10, 55];
    min_n_pack_feat_pnl = 5; 
end

f = figure('color','white', ...
           'position', fig_size); 
       
pnl = panel(f); 

pnl.pack('v', [20, 80]); 

pnl(2).pack('v', n_cond); 


ymax_mX = -Inf; 
ymax_mX_subtr = -Inf; 

colnames = tbl.Properties.VariableNames; 
datanames = fieldnames(data_to_plot); 

for i_cond=1:n_cond

    fprintf('plotting %d/%d\n', i_cond, n_cond);

    x = data_to_plot(i_cond).x; 
    t = data_to_plot(i_cond).t; 
    mX = data_to_plot(i_cond).mX; 
    freq = data_to_plot(i_cond).freq; 
    acf = data_to_plot(i_cond).acf; 
    lags = data_to_plot(i_cond).lags; 

    ap = []; 
    mX_subtr = []; 
    acf_subtr = []; 

    if plot_subtr
        if any(strcmp(datanames, 'ap'))
            ap = data_to_plot(i_cond).ap; 
        end
        if any(strcmp(datanames, 'mX_subtr'))
            mX_subtr = data_to_plot(i_cond).mX_subtr; 
        end
        if any(strcmp(datanames, 'acf_subtr'))
            acf_subtr = data_to_plot(i_cond).acf_subtr; 
        end
    end
    
    % update yaxis maximum for FFT
    frex_idx = dsearchn(freq', par.frex');
    ymax_mX = max(ymax_mX, max(mX(frex_idx)));
    
    if ~isempty(mX_subtr)
        ymax_mX_subtr = max(ymax_mX_subtr, max(mX_subtr(frex_idx)));
    end
    
    plot_example(x, t, acf, lags, ap, mX, freq, ...
                 par.lags_meter_rel, par.lags_meter_unrel, ...
                 par.freq_meter_rel, par.freq_meter_unrel, ...
                 'mX_subtr', mX_subtr, ...
                 'acf_subtr', acf_subtr, ...
                 'pnl', pnl(2, i_cond), ...
                 'subplot_proportions', example_subplot_proportions, ...
                 'max_lag', par.max_lag, ...
                 'max_freq', par.max_freq_plot, ...
                 'plot_time_xaxis', i_cond == n_cond, ...
                 'plot_xlabels', i_cond == n_cond, ...
                 'plot_xticks', i_cond == n_cond, ...
                 'plot_features', false, ...
                 'time_col', colors{i_cond}, ...
                 'prec', 1e6, ...
                 'fontsize', par.fontsize, ...
                 'normalize_acf_for_plotting', false, ...
                 varargin_for_examples{:});              
             
    if ~isempty(mX_subtr) || ~isempty(acf_subtr)
        pnl(2, i_cond).margintop = 25; 
    else
        pnl(2, i_cond).margintop = 15; 
    end
end


for i_cond=1:n_cond

    ax = pnl(2, i_cond, 2, 1).select();
    ax.YLim = [0, ymax_mX];
    
    if ~isempty(mX_subtr)
        ax = pnl(2, i_cond, 2, 2).select();
        ax.YLim = [0, ymax_mX_subtr];
    end
    
end

%%

prec = 1e2; 

ytick_at_means = false;

pnl(1).pack('h', max(min_n_pack_feat_pnl, length(feat_to_plot))); 

for i_feat=1:length(feat_to_plot)
    
    switch feat_to_plot{i_feat}
        
        case 'z_meter_acf'
            feat_label = 'zscore'; 
            tit = 'ACF'; 
            feat_orig = []; 
            
        case 'z_meter_fft'
            feat_label = 'zscore'; 
            tit = 'FFT'; 
            feat_orig = []; 
            
        case 'z_meter_acf_raw'
            feat_label = 'zscore (raw)'; 
            tit = 'ACF'; 
            feat_orig = 'z_meter_acf_orig'; 
            
        case 'z_meter_fft_raw'
            feat_label = 'zscore (raw)'; 
            tit = 'FFT'; 
            feat_orig = 'z_meter_fft_orig'; 

        case 'z_meter_acf_subtr'
            feat_label = 'zscore (1/f-corrected)'; 
            tit = 'ACF'; 
            feat_orig = 'z_meter_acf_orig'; 
            
        case 'z_meter_fft_subtr'
            feat_label = 'zscore (1/f-corrected)'; 
            tit = 'FFT'; 
            feat_orig = 'z_meter_fft_orig'; 
            
        case 'ap_offset'
            feat_label = 'offset'; 
            tit = 'AP'; 
            feat_orig = [];             
            
        case 'ap_exponent'
            feat_label = 'exponent'; 
            tit = 'AP'; 
            feat_orig = [];             
    end
    
    % check 
    if ~any(strcmp(colnames, feat_orig))
       warning('column with ground-truth (orig) feature not present...skipping');  
       feat_orig = []; 
    end
    
    pnl(1, i_feat).pack({[0, 0 , 1, 1]}); 
    ax = pnl(1, i_feat, 1).select(); 
    
    plot_multiple_cond(tbl, ...
                      'cond_column', cond_colname, ...
                      'feat_column', feat_to_plot{i_feat}, ...
                      'feat_orig_column', feat_orig, ...
                      'ax', ax, ...
                      'colors', colors, ...
                      'plot_legend', true, ...
                      varargin_for_points{:}); 
                      
    ax.YLim = [min(min(floor(tbl{:, feat_to_plot{i_feat}}*prec)/prec), Inf), ...
               max(max(ceil(tbl{:, feat_to_plot{i_feat}}*prec)/prec), -Inf), ...
               ];
    
    ax.YTick = ax.YLim; 
       
    pnl(1, i_feat).ylabel(sprintf('%s', feat_label)); 
    
    % make the figure nice 
    pnl(1, i_feat).margin = [19, 5, 5, 40]; 

    pnl(1, i_feat).title(tit); 

    pnl(1, i_feat).fontsize = par.fontsize; 

end

% fix legend position
sw = 1;
for i=1:length(f.Children)
    if strcmpi(f.Children(i).Type, 'legend') 
       if sw
           f.Children(i).Title.Interpreter = 'none'; 
           f.Children(i).Title.String = cond_colname; 
           f.Children(i).Position(1) = 0.94; 
           f.Children(i).Position(2) = 0.85; 
           sw = 0;
       else
           f.Children(i).Visible = 'off';
       end
    end
end

pnl(1).de.marginright = 40;
pnl(2).margintop = 15;
pnl.margin = [25, 10, 25, 15];












