clear

par = get_par(); 

cmap_name = 'RdPu'; 

cond_colname = 'duty_cycle'; 

feat_to_plot = {
    'z_meter_fft'
    'z_meter_acf'
    }; 

fname = '03_ir';

f = plot_multi_figure(fname, cmap_name, cond_colname, feat_to_plot);

save_fig(f, fullfile(par.data_path, [fname, '.svg']));
