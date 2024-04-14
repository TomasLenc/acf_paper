% Plots to demontrate sensitivity of ACF to shift and offset of the input
% signal (see Figure 6 in the paper)

clear

par = get_par(); 

% make clean signal for the whole trial 
[x_clean, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', par.ir ...
                    );
                
noise = get_colored_noise2(...
    [1, round(par.trial_dur*par.fs)], ...
    par.fs, -1.5); 

x = add_signal_noise(x_clean, noise, 2);

x = zscore(x); 

lags_meter_rel = [0.8]; 

lags_meter_unrel = [0.4, 0.6, 1.0, 1.2]; 

xlims = [0, 1.6]; 

%%

n = 10; 

% ==================================
param_name = 'scale'; % scale, shift
% ==================================

if strcmp(param_name, 'scale')
    
    param_vals = linspace(1, 3, n); 
    % param_vals = logspace(log10(1), log10(n), n); 
    fun = @(x, b) b * x;  
    cmap_name = 'Greys'; 
    
elseif strcmp(param_name, 'shift')
    
    param_vals = [1:n] / n * range(x)/2; 
    % param_vals = logspace(log10(1/n), log10(1), n); 
    fun = @(x, a) a + x;  
    cmap_name = 'Blues'; 
    
end

% =========================================================================

colors = num2cell(brewermap(n + round(n/3), cmap_name), 2); 
colors = colors(end-n+1:end); 

f = figure('color','white', 'position', [356 219 716 546]); 
pnl = panel(f); 
pnl.pack('h', [80, 20])
pnl(1).pack({[0, 0, 1, 1]}); 

ax = pnl(1, 1).select(); 

r_mean = nan(1, length(param_vals)); 
r_ratio = nan(1, length(param_vals)); 
acf_mean = nan(1, length(param_vals)); 
acf_ratio = nan(1, length(param_vals)); 
acf_contrast = nan(1, length(param_vals)); 
acf_z = nan(1, length(param_vals)); 

for i_cond=1:length(param_vals)
    
    x_current = fun(x, param_vals(i_cond)); 

    [acf, lags] = get_acf(x_current, par.fs, ...
                       'normalize_x', false, ...
                       'force_x_positive', false, ...
                       'normalize_acf_to_1', false, ...
                       'normalize_acf_z', false ...
                       );    
                   
    plot(lags, acf, 'linew', 1.4, 'color', colors{i_cond}); 
    hold on
    
    % Notice that when we just take Pearson correlation with a lagged version
    % of the signal, this is invariant to shit and scale 
    lags_meter_rel_idx = dsearchn(lags', lags_meter_rel')'; 
    lags_meter_unrel_idx = dsearchn(lags', lags_meter_unrel')'; 

    r_meter_rel = nan(1, length(lags_meter_rel));
    for i=1:length(lags_meter_rel)
        r = corrcoef(x_current, circshift(x_current, lags_meter_rel_idx(i))); 
        r_meter_rel(i) = r(2); 
    end
    r_meter_unrel = nan(1, length(lags_meter_unrel));
    for i=1:length(lags_meter_unrel)
        r = corrcoef(x_current, circshift(x_current, lags_meter_unrel_idx(i))); 
        r_meter_unrel(i) = r(2); 
    end
    
    
    r_mean(i_cond) = mean(r_meter_rel);  
        
    r_ratio(i_cond) = mean(r_meter_rel) / mean(r_meter_unrel); 
    
    acf_meter_rel = acf(lags_meter_rel_idx); 
    acf_meter_unrel = acf(lags_meter_unrel_idx); 

    acf_mean(i_cond) = mean(acf_meter_rel); 

    acf_ratio(i_cond) = mean(acf_meter_rel) / mean(acf_meter_unrel); 

    acf_contrast(i_cond) = (mean(acf_meter_rel) - mean(acf_meter_unrel)) /...
                           (mean(acf_meter_rel) + mean(acf_meter_unrel)); 
    
    z = zscore([acf_meter_rel, acf_meter_unrel]);
    acf_z(i_cond) = mean(z(1:length(acf_meter_rel))); 
end


%% 

prec = 1e4; 

ylims = ax.YLim; 
y_to_plot = [floor(ylims(1)*prec)/prec, ceil(ylims(2)*prec)/prec]; 
h = plot([lags_meter_rel; lags_meter_rel], y_to_plot, ...
         'linew', 3, 'color', [0.8706    0.1765    0.1490]); 
for i=1:length(h)
     h(i).Color(4) = 0.3; 
end
h = plot([lags_meter_unrel; lags_meter_unrel], y_to_plot, ...
         'linew', 3, 'color',  [0.1922    0.5098    0.7412]); 
for i=1:length(h)
     h(i).Color(4) = 0.3; 
end

ax.YLim = ylims; 
ax.XLim = xlims; 

cbar = colorbar('WestOutside'); 
colormap(cell2mat(colors)); 
cbar.Ticks = linspace(0+0.5/length(param_vals), 1-0.5/length(param_vals), length(param_vals)); 
cbar.TickLabels = cellfun(@(x) round(x, 2), num2cell(param_vals)); 
cbar.Label.String = param_name; 
cbar.Position(1) = cbar.Position(1) - abs(cbar.Position(1)) * 0.2; 
cbar.Position(2) = cbar.Position(2) +  abs(cbar.Position(2)) * 1; 
cbar.Position(4) = cbar.Position(4) * 0.75; 

ax.XAxisLocation = 'origin'; 
ax.YTick = ax.YTick([1, end]);
ax.XTick = []; 
ax.TickDir = 'out'; 

pnl(1).ylabel('autocorrelation')




pnl(2).pack('v', 5); 

prec = 1e9; 
ylim_buffer = 0.5; 


ax = pnl(2, 1).select(); 
scatter(param_vals, acf_mean, 15, cell2mat(colors), 'filled')
ylims = [floor(ax.YLim(1)*prec)/prec, ...
         ceil(ax.YLim(2)*prec)/prec]; 
ylims = [ylims(1) - range(ylims)*ylim_buffer, ylims(2) + range(ylims)*ylim_buffer]; 
ax.YLim = ylims; 
ax.XTick = []; 
ax.YTick = []; 
ax.Title.String = 'raw acf'; 

ax = pnl(2, 2).select(); 
scatter(param_vals, r_mean, 15, cell2mat(colors), 'filled')
ylims = [floor(ax.YLim(1)*prec)/prec, ...
         ceil(ax.YLim(2)*prec)/prec]; 
ylims = [ylims(1) - range(ylims)*ylim_buffer, ylims(2) + range(ylims)*ylim_buffer]; 
ax.YLim = ylims; 
ax.XTick = []; 
ax.YTick = []; 
ax.Title.String = 'Pearson r'; 

ax = pnl(2, 3).select(); 
scatter(param_vals, acf_ratio, 15, cell2mat(colors), 'filled')
ylims = [floor(ax.YLim(1)*prec)/prec, ...
         ceil(ax.YLim(2)*prec)/prec]; 
ylims = [ylims(1) - range(ylims)*ylim_buffer, ylims(2) + range(ylims)*ylim_buffer]; 
ax.YLim = ylims; 
ax.XTick = []; 
ax.YTick = []; 
ax.Title.String = 'ratio of acf'; 

ax = pnl(2, 4).select(); 
scatter(param_vals, acf_contrast, 15, cell2mat(colors), 'filled')
ylims = [floor(ax.YLim(1)*prec)/prec, ...
         ceil(ax.YLim(2)*prec)/prec]; 
ylims = [ylims(1) - range(ylims)*ylim_buffer, ylims(2) + range(ylims)*ylim_buffer]; 
ax.YLim = ylims; 
ax.XTick = []; 
ax.YTick = []; 
ax.Title.String = 'contrast of acf'; 

ax = pnl(2, 5).select(); 
scatter(param_vals, acf_z, 15, cell2mat(colors), 'filled')
ylims = [floor(ax.YLim(1)*prec)/prec, ...
         ceil(ax.YLim(2)*prec)/prec]; 
ylims = [ylims(1) - range(ylims)*ylim_buffer, ylims(2) + range(ylims)*ylim_buffer]; 
ax.YLim = ylims; 
ax.XTick = []; 
ax.YTick = []; 
ax.Title.String = 'zscore of acf'; 

pnl(2).marginleft = 10; 
pnl.margin = [45, 10, 5, 7]; 
pnl.fontsize = 16; 

%% 

save_path = fullfile(par.fig_path, 'general', 'acf_normalization'); 
mkdir(save_path); 

fname = sprintf('acf_normalization_%s', param_name); 
print(fullfile(save_path, fname), '-dsvg', '-painters', f);  
print(fullfile(save_path, fname), '-dpng', '-painters', '-r300', f);  

%% 

% print resutls 
fprintf('\n------------------------ %s ----------------------\n', param_name);
disp(struct2table(struct(param_name, param_vals', ...
                         'Pearson_r', r_mean', ...
                         'ratio_of_Pearson_r', r_ratio', ...
                         'ratio_of_ACF', acf_ratio',...
                         'zscore_of_ACF', acf_z')))

                



%%

f = figure('color','white', 'position',[356 599 835 166]); 
pnl = panel(f);
pnl.pack('h', 2); 
        

ax = pnl(1).select(); 

x_current = fun(x, param_vals(1)); 
plot(ax, t, x_current, 'color', colors{1}, 'linew', 2); 
hold(ax, 'on');  

ax.XLim = [0, 4.8]; 
ax.XAxis.Visible = 'off'; 
ax.YAxis.Visible = 'off'; 
ylims = ax.YLim; 


ax = pnl(2).select(); 

x_current = fun(x, param_vals(end)); 
plot(ax, t, x_current, 'color', colors{end}, 'linew', 2); 
hold(ax, 'on');  

ax.XLim = [0, 4.8]; 
ax.XAxis.Visible = 'off'; 
ax.YAxis.Visible = 'off'; 
ylims(1) = min(ylims(1), ax.YLim(1)); 
ylims(2) = max(ylims(2), ax.YLim(2)); 


ax = pnl(1).select(); 
ax.YLim = ylims; 
ax.YTick = []; 

ax = pnl(2).select(); 
ax.YLim = ylims; 
ax.YTick = []; 

pnl.de.margin = 1; 
pnl.margin = 5; 

fname = sprintf('acf_normalization_%s_timeDomain', param_name); 
print(fullfile(save_path, fname), '-dsvg', '-painters', f);  
print(fullfile(save_path, fname), '-dpng', '-painters', '-r300', f);  


                               