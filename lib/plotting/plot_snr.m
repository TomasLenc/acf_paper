function plot_snr(par)
  
fname = sprintf('ir-%s_noise-%s_apFitMethod-%s_onlyHarm-%s_snr', ...
               par.ir_type, ...
               par.noise_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 
           
cmap_name = 'OrRd'; 

cond_colname = 'snr'; 

feat_to_plot = {
    'z_meter_fft_raw'
    'z_meter_fft_subtr'
    'z_meter_acf_raw'
    'z_meter_acf_subtr'
    }; 

[f, pnl] = plot_multi_figure(par.data_path, fname, cmap_name, cond_colname, feat_to_plot, ...
                             'varargin_for_points', {'zero_line', true});

all_ax = {pnl(1,1,1).select(), pnl(1,2,1).select()}; 
match_ylims(all_ax); 

all_ax = {pnl(1,3,1).select(), pnl(1,4,1).select()}; 
match_ylims(all_ax); 

save_fig(f, fullfile(par.data_path, [fname, '.svg']));


%%
