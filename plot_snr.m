clear

par = get_par(); 

cmap_name = 'OrRd'; 

cond_colname = 'snr'; 

feat_to_plot = {
    'z_meter_fft_raw'
    'z_meter_fft_subtr'
    'z_meter_acf_raw'
    'z_meter_acf_subtr'
    }; 

fname = '04_snr';

[f, pnl] = plot_multi_figure(fname, cmap_name, cond_colname, feat_to_plot);

all_ax = {pnl(1,1,1).select(), pnl(1,2,1).select()}; 
match_ylims(all_ax); 

all_ax = {pnl(1,3,1).select(), pnl(1,4,1).select()}; 
match_ylims(all_ax); 

save_fig(f, fullfile(par.data_path, [fname, '.svg']));


%%
