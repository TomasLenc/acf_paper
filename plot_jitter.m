function plot_jitter(par)
  
fname = sprintf('ir-%s_apFitMethod-%s_onlyHarm-%s_jitter', ...
               par.ir_type, ...
               par.ap_fit_method, ...
               jsonencode(par.only_use_f0_harmonics)); 
           
cmap_name = 'YlGn'; 

cond_colname = 'jitter'; 

feat_to_plot = {
    'z_meter_fft_raw'
    'z_meter_acf_raw'
    }; 

[f, pnl] = plot_multi_figure(par.data_path, fname, cmap_name, cond_colname, feat_to_plot, ...
                             'plot_subtr', false, ...
                             'varargin_for_points', {'zero_line', true});
 
for i=1:length(pnl(2).children)
    pnl(2, i).margintop = 6; 
end
% f.Position = [334 1202 1319 974]; 


save_fig(f, fullfile(par.data_path, [fname, '.svg']));


%%
