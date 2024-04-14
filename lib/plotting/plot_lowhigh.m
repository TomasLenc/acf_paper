function plot_lowhigh(par, varargin)

do_chunk_acf = false; 
if any(strcmp(varargin, 'chunk_acf'))
    do_chunk_acf = true; 
end

do_chunk_time = false; 
if any(strcmp(varargin, 'chunk_time'))
    do_chunk_time = true; 
end


if do_chunk_time
    fname = sprintf('exp-lowhigh_apFitMethod-%s_onlyHarm-%s_roi-%s_chunk-*_eegIndividual.mat', ...
                    par.ap_fit_method, ...
                    jsonencode(par.only_use_f0_harmonics),...
                    par.roi_name); 
else
    fname = sprintf('exp-lowhigh_apFitMethod-%s_onlyHarm-%s_roi-%s_eegIndividual.mat', ...
                    par.ap_fit_method, ...
                    jsonencode(par.only_use_f0_harmonics),...
                    par.roi_name); 
end
d = dir(fullfile(par.data_path, fname)); 
[~, fname] = fileparts(d.name);     

tmp = load(fullfile(par.data_path, [fname, '.mat'])); 
data_to_plot = tmp.data_to_plot; 

rhythms = {'unsyncopated', 'syncopated'}; 
tones = {'L', 'H'}; 

% colors = {
%     [181, 44, 13]/255 
%     [16, 88, 176]/255
%     }; 

colors = {
    [0 0 0 ]/255 
    [0 0 0]/255
    }; 

%%

fig_pos = [155 335 1657 484]; 
subplot_ratio = [15, 85]; 

if do_chunk_acf || do_chunk_time || par.max_lag < 5
    fig_pos = [249 201 600 450];
    subplot_ratio = [40, 60]; 
end


f = figure('color', 'white', 'Position', fig_pos); 
pnl = panel(f); 

pnl.pack('v', 2); 

for i_tone=1:2
    pnl(i_tone).pack('h', 2); 
    for i_rhythm=1:2
        pnl(i_tone, i_rhythm).pack('v', 2); 
        
        pnl(i_tone, i_rhythm, 1).pack('h', subplot_ratio); 
        pnl(i_tone, i_rhythm, 2).pack('h', subplot_ratio); 
    end
end
% pnl.select('all'); 

pnl.de.margin = [5, 1, 3, 5]; 
pnl(2).margintop = 15; 
pnl(1, 2).marginleft = 15; 
pnl(2, 2).marginleft = 15; 
pnl.margin = [15, 10, 5, 10]; 

xlim_acf = [0, min(par.max_lag, 13.2)]; 
xtick_acf = [0, min(par.max_lag, 13.2)]; 
linew_acf = 1; 
prec = 100; 

    
for i_tone=1:2

    tone = tones{i_tone}; 

    pnl(i_tone).title(tone); 

    for i_rhythm=1:2

        rhythm = rhythms{i_rhythm}; 

        pnl(i_tone, i_rhythm).ylabel(rhythm); 
    
                
        mask = strcmp({data_to_plot.rhythm}, rhythm) & ...
               strcmp({data_to_plot.tone}, tone) ; 
        
           
        freq_coch = data_to_plot(mask).freq_coch; 
        mX_coch = data_to_plot(mask).mX_coch; 
        freq = data_to_plot(mask).freq; 
        mX_subtr = data_to_plot(mask).mX_subtr; 
           
        lags = data_to_plot(mask).lags; 
        acf_subtr = data_to_plot(mask).acf_subtr; 
        
        lags_coch = data_to_plot(mask).lags_coch; 
        acf_coch = data_to_plot(mask).acf_coch; 
        
        lags_meter_rel = par.lags_meter_rel; 
        lags_meter_unrel = par.lags_meter_unrel; 
        
        % average acf cycles 
        if do_chunk_acf
            acf_coch = epoch_chunks(acf_coch, 1/(lags_coch(2) - lags_coch(1)), 2.4); 
            acf_coch = mean(acf_coch, 1); 
            lags_coch = lags_coch(1 : length(acf_coch)); 

            acf_subtr = epoch_chunks(acf_subtr, 1/(lags(2) - lags(1)), 2.4); 
            acf_subtr = squeeze(mean(acf_subtr, 1)); 
            lags = lags(1 : length(acf_subtr)); 
            
            lags_meter_rel = [0.8, 1.6]; 
            lags_meter_unrel = []; 
        end
        
        % coch
        ax = pnl(i_tone, i_rhythm, 1, 1).select(); 
        
        plot_fft(freq_coch, mX_coch, ...
         'ax', ax, ...
         'frex_meter_rel', par.freq_meter_rel, ...
         'frex_meter_unrel', par.freq_meter_unrel, ...
         'maxfreqlim', par.max_freq); 
        
        ax.YTick = []; 
        
        
        ax = pnl(i_tone, i_rhythm, 1, 2).select(); 
        
        acf_to_plot = zscore(acf_coch, [], 2); 
        acf_to_plot = acf_to_plot + min(acf_to_plot); 

        plot_acf(ax, ...
                 acf_to_plot, ...
                 lags_coch, ...
                 'lags_meter_rel', lags_meter_rel, ...
                 'lags_meter_unrel', lags_meter_unrel, ...
                 'linew', linew_acf, ...
                 'opacity_lagz', 0.5, ...
                 'prec', prec); 
        ax.XLim = xlim_acf;  
        ax.YTick = []; 
        ax.XAxis.Visible = 'off'; 

        
        % EEG 
        ax = pnl(i_tone, i_rhythm, 2, 1).select(); 
        
        plot_fft(freq, mean(mX_subtr, 1), ...
         'ax', ax, ...
         'frex_meter_rel', par.freq_meter_rel, ...
         'frex_meter_unrel', par.freq_meter_unrel, ...
         'maxfreqlim', par.max_freq); 
                
        
        ax = pnl(i_tone, i_rhythm, 2, 2).select(); 

        N = size(acf_subtr, 1); 
        acf_subtr_norm = zscore(acf_subtr, [], 2); 
        
        acf_to_plot = mean(acf_subtr_norm, 1); 
        offset = min(acf_to_plot); 
        acf_to_plot = acf_to_plot + offset; 

        sem_to_plot = std(acf_subtr_norm, [], 1) / sqrt(N) + offset; 
        ci_to_plot = sem_to_plot * norminv(1-0.05); 
        
        plot_acf(ax, ...
                 acf_to_plot, ...
                 lags, ...
                 'lags_meter_rel', lags_meter_rel, ...
                 'lags_meter_unrel', lags_meter_unrel, ...
                 'linew', linew_acf, ...
                 'col_acf', colors{i_tone}, ...
                 'opacity_lagz', 0.5, ...
                 'prec', prec); 
        ax.XLim = xlim_acf;  
        ax.YTick = []; 
        ax.XAxis.Visible = 'off'; 
        
        hold(ax, 'on'); 
        
        error_area = [acf_to_plot - sem_to_plot, flip(acf_to_plot + sem_to_plot)]; 
        
        fill(ax, [lags, flip(lags)], ...
             error_area, ...
             colors{i_tone}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
         
         ax.YLim = [min(error_area), max(error_area)]; 
        
    end
end

ax.XAxis.Visible = 'on'; 
ax.XTick = xtick_acf; 
prec = 10; 
ax.XTickLabel = xtick_acf; 

pnl.fontsize = 12; 


% fix ylimx for FFT
ymax_fft = -Inf; 
for i_tone=1:2
    for i_rhythm=1:2
        ax = pnl(i_tone, i_rhythm, 2, 1).select(); 
        ymax_fft = max(ymax_fft, max([ax.Children(1:2).YData])); 
    end
end
for i_tone=1:2
    for i_rhythm=1:2
        ax = pnl(i_tone, i_rhythm, 2, 1).select(); 
        ax.YLim = [0, ymax_fft]; 
        ax.YTick = [0, floor(ymax_fft * 100) / 100]; 
    end
end



%%

fname_to_save = sprintf('%s_fft-acf', fname); 

if do_chunk_acf 
    fname_to_save = sprintf('%s_chunk_acf', fname); 
end

save_fig(f, fullfile(par.data_path, fname_to_save));



