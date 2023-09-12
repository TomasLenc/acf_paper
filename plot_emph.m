
clear

par = get_par(); 

cmap_name = 'YlGn'; 

cond_colname = 'emphasis'; 

feat_to_plot = {
    'z_meter_fft'
    'z_meter_acf'
    }; 

fname = '05_emph';

f = plot_multi_figure(fname, cmap_name, cond_colname, feat_to_plot);

save_fig(f, fullfile(par.data_path, [fname, '.svg']));
