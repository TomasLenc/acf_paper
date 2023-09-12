
clear

par = get_par(); 

cmap_name = 'YlGn'; 

cond_colname = 'jitter'; 

feat_to_plot = {
    'z_meter_fft_raw'
    'z_meter_acf_raw'
    }; 

fname = '10_jitter';

[f, pnl] = plot_multi_figure(fname, cmap_name, cond_colname, feat_to_plot, ...
                             'plot_subtr', false, ...
                             'varargin_for_points', {'zero_line', true});
 
for i=1:length(pnl(2).children)
    pnl(2, i).margintop = 6; 
end
% f.Position = [334 1202 1319 974]; 


save_fig(f, fullfile(par.data_path, [fname, '.svg']));


%%
