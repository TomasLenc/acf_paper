function main03_ir(par, varargin)
% clear 
% varargin = {}; 
% par = get_par(); 

parser = inputParser; 

addParameter(parser, 'ir_type', 'square'); % square, erp, erp2

parse(parser, varargin{:});

ir_type = parser.Results.ir_type;

addpath(genpath(par.acf_tools_path)); 
addpath(genpath(par.rnb_tools_path)); 
addpath(genpath('lib'))


%% simulate

noise_exponent = -1.5; 

fit_knee = false; 

% percent extreme values omitted for plotting
ylim_quantile_cutoff = 0.05; 

% ------------------------------------------------
if strcmp(ir_type, 'square')
    cond_type = 'duty cycle';
    duty_cycles = linspace(0.050, 0.180, 6); 
end
if strcmp(ir_type, 'erp')
    cond_type = 'IR freq'; 
    duty_cycles = linspace(10, 4, 6); 
end
% ------------------------------------------------


%%

n_cond = length(duty_cycles); 

% colors
cmap_name = '-RdPu'; 
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(1:n_cond, :); 


%% test performance across duty cycles


% allocate 
feat_acf_orig = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
    'contrast_meter_rel', []); 

feat_acf = struct(...
    'z_meter_rel', [], 'ratio_meter_rel', [], ...
    'ratio_meter_rel_left', [], 'ratio_meter_rel_right', [],...
    'contrast_meter_rel', []); 


feat_fft_orig = struct('z_meter_rel', []); 

feat_fft = struct('z_meter_rel', []); 

cond_labels = {}; 

ymax_mX = -Inf;

%%

f = figure('color','white', ...
           'position', [146 1340 1329 933]); 
       
pnl = panel(f); 

pnl.pack('v', [20, 80]); 

pnl(2).pack('v', n_cond); 

example_subplot_proportions = [35, 10, 55];



for i_cond=1:n_cond
    
    duty_cycle = duty_cycles(i_cond); 
        
    cond_labels{i_cond} = sprintf('%.2f', duty_cycle); 

    fprintf('calculating %d/%d\n', i_cond, n_cond)

    if strcmp(ir_type, 'erp')
        ir = get_erp_kernel(par.fs,...
            'amplitudes', 1,...
            't0s', 0, ...
            'taus', 0.050, ...
            'f0s', duty_cycle, ...
            'duration', 0.2 ...
            ); 
    elseif strcmp(ir_type, 'square')
        ir = get_square_kernel(par.fs, ...
            'duration', duty_cycle, ...
            'rampon', 0, ...
            'rampoff', 0 ...
            ); 
    else
        error('ir kind not recognized'); 
    end
    
    % make whole signal 
    [x, t] = get_s(...
                        par.pat, ...
                        par.grid_ioi, ...
                        par.fs, ...
                        'n_cycles', par.n_cycles, ...
                        'ir', ir ...
                        );
    
    % get acf
    % -------
                                       
    % withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, par.fs);    
                                                                          
    % get features
    % ------------
    feat_acf(i_cond) = get_acf_features(acf, lags, ...
                                par.lags_meter_rel, par.lags_meter_unrel);    
                                                                  
    feat_fft(i_cond) = get_fft_features(mX, freq, ...
                                            par.freq_meter_rel, par.freq_meter_unrel); 
                             
    % plot example 
    % ------------
    rep_to_plot_idx = 1; 
    
    % update yaxis maximum for FFT
    frex_idx = dsearchn(freq', par.frex');
    amps = mX(rep_to_plot_idx, frex_idx);
    ymax_mX = max(ymax_mX, max(amps));

    plot_example(x(rep_to_plot_idx, :), t, ...
                     acf(rep_to_plot_idx, :), lags, ...
                     [], ...
                     mX(rep_to_plot_idx, :), freq, ...
                     par.lags_meter_rel, par.lags_meter_unrel, ...
                     par.freq_meter_rel, par.freq_meter_unrel, ...
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
                     'normalize_acf_for_plotting', false);                                        
    f.Name = cond_labels{i_cond};     
    pnl(2, i_cond).margintop = 15; 

end


for i_cond=1:n_cond
    ax = pnl(2, i_cond, 2, 1).select();
    ax.YLim = [0, ymax_mX];
end

% assign labels
[feat_acf.name] = deal(cond_labels{:}); 
[feat_fft.name] = deal(cond_labels{:}); 

% assign colors
[feat_acf.color] = deal(colors{:}); 
[feat_fft.color] = deal(colors{:}); 


%%


% plot [35, 10, 55]
% ----
cond_to_plot = {
    'acf-z_meter_rel'
    'fft-z_meter_rel'
    }; 

pnl(1).pack('h', 5); 

for i_cond=1:length(cond_to_plot)
    
    ytick_at_means = false;
    
    switch cond_to_plot{i_cond}
        
        case 'acf-mean_meter_rel'
            feat_raw = feat_acf; 
            feat_fieldname = 'mean_meter_rel'; 
            feat_label = 'mean'; 
            tit = 'ACF'; 
        case 'acf-ratio_meter_rel'
            feat_raw = feat_acf; 
            feat_fieldname = 'ratio_meter_rel'; 
            feat_label = 'ratio'; 
            tit = 'ACF'; 
        case 'acf-z_meter_rel'
            feat_raw = feat_acf; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'ACF'; 
        case 'fft-z_meter_rel'
            feat_raw = feat_fft; 
            feat_fieldname = 'z_meter_rel'; 
            feat_label = 'zscore'; 
            tit = 'FFT'; 

    end
    
    
    
    pnl(1, i_cond).pack({[0, 0 , 1, 1]}); 
    ax = pnl(1, i_cond, 1).select(); 

    feat = RenameField(feat_raw, feat_fieldname, 'data');

    plot_multiple_cond('ax', ax, ...
                      'plot_legend', true, ...
                      'ytick_at_means', ytick_at_means, ...
                      'feat', feat, ...
                      'ylim_quantile_cutoff', ylim_quantile_cutoff); 

    pnl(1, i_cond).ylabel(sprintf('%s raw', feat_label)); 

    
    ax.YLim = [min(-max(abs([feat.data])), -0.1), ...
               max(max(abs([feat.data])), 0.1), ...
               ];
    
    ax.YTick = ax.YLim; 
    
    % make the figure nice 
    pnl(1, i_cond).margin = [19, 5, 5, 40]; 

    pnl(1, i_cond).title(tit); 

    pnl(1, i_cond).fontsize = par.fontsize; 


    

end

% fix legend position
sw = 1;
for i=1:length(f.Children)
    if strcmpi(f.Children(i).Type, 'legend') 
       if sw
           f.Children(i).Title.String = cond_type; 
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

fname = sprintf('03_ir_irType-%s_exp-%.1f_%s_%s', ...
                 ir_type, noise_exponent, tit, feat_label);
if par.save_figs
   save_fig(f, fullfile(par.fig_path, fname))
end

% save parameters 
save(fullfile(par.fig_path, [fname, '_par.mat']), 'par'); 