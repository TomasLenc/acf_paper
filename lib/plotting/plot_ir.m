function plot_ir(par)

fname = sprintf('ir-%s_ir', par.ir_type);

cmap_name = 'RdPu'; 

cond_colname = 'duty_cycle'; 

feat_to_plot = {
    'z_meter_fft'
    'z_meter_acf'
    }; 

f = plot_multi_figure(par.data_path, fname, cmap_name, cond_colname, feat_to_plot);

save_fig(f, fullfile(par.data_path, [fname, '.svg']));
