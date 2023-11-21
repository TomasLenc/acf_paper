clear

par = get_par(); 

%%

n_pats = 5; 

n_events = 12; 
n_sounds = 8; 
max_group_size = 4; 

all_good_pats = find_all_patterns(n_events, n_sounds, max_group_size); 

% idx = randsample(size(all_good_pats, 1), n_pats); 

par.all_pats = all_good_pats(:, :); 

%%


f = figure('color','white', 'position', [674 33 1000 923]); 
pnl = panel(f); 
pnl.pack('v', size(par.all_pats, 1)+1);

for i=1:size(par.all_pats, 1) + 1
    pnl(i).pack('h', 6);
end

pnl.de.margin = [5, 1, 5, 1]; 

for i_pat=1:size(par.all_pats, 1)

    pat = repmat(par.all_pats(i_pat, :), 1, 1); 

    
    % ========================================================================
    
    ax = pnl(i_pat, 1).select(); 
    t_pattern = [0 : length(pat)-1] * par.grid_ioi; 
    scatter(t_pattern(logical(pat)), zeros(1, sum(pat)), 70, 'kx', 'LineWidth', 2)
    hold on
    scatter(t_pattern(logical(~pat)), zeros(1, sum(~pat)), 30, 'k.', 'LineWidth', 2)
    ax.XLim = [0, length(pat) * par.grid_ioi]; 
    ax.Visible = 'off'; 
    
    
    ax = pnl(i_pat, 2).select(); 

    mX = abs(fft(pat))  / length(pat); 
    bar(ax, mX)
    ax.XTickLabel = [0 : length(pat)-1]; 
    ax.YAxis.Visible = 'off';
    
    % ========================================================================
    
    
    ir = 1;
    
    N = round(par.grid_ioi * length(pat) * par.fs); 
    x = zeros(1, N); 
    for i=1:length(pat)
        if pat(i)
            idx = round((i-1) * par.grid_ioi * par.fs); 
            x(idx+1:idx+length(ir)) = ir; 
        end
    end

    ax = pnl(i_pat, 3).select(); 
    plot(ax, x, 'linew', 2, 'color', 'k'); 
    ax.Visible = 'off'; 
    ax.YLim = [-1, 2];
    
    

    ax = pnl(i_pat, 4).select(); 
    mX = abs(fft(x))  / length(x); 
    freq = [0 : N-1] / N * par. fs; 
    plot_fft(freq, mX, 'maxfreqlim', 15, 'linew', 2, 'ax', ax);
    ax.XTick = [0, 5, 10, 15]; 
    if i_pat == size(par.all_pats, 1)
        ax.XAxis.Visible = 'on';
    end
    ax.YAxis.Visible = 'off';
        
    
    % ========================================================================
    
    
    ir = get_square_kernel(par.fs, ...
            'duration', 0.100, ...
            'rampon', 0, ...
            'rampoff', 0 ...
            ); 
    
    N = round(par.grid_ioi * length(pat) * par.fs); 
    x = zeros(1, N); 
    for i=1:length(pat);
        if pat(i)
            idx = round((i-1) * par.grid_ioi * par.fs); 
            x(idx+1:idx+length(ir)) = ir; 
        end
    end

    ax = pnl(i_pat, 5).select(); 
    plot(ax, x, 'linew', 2, 'color', 'k'); 
    ax.Visible = 'off'; 
    ax.YLim = [-1, 2];

    
    
    ax = pnl(i_pat, 6).select(); 
    mX = abs(fft(x))  / N; 
    freq = [0 : N-1] / N * par. fs; 
    plot_fft(freq, mX, 'maxfreqlim', 15, 'linew', 2, 'ax', ax);
    ax.XTick = [0, 5, 10, 15]; 
    if i_pat == size(par.all_pats, 1)
        ax.XAxis.Visible = 'on';
    end
    ax.YAxis.Visible = 'off';
    
    
end

ir = get_square_kernel(par.fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 

ir = [ir, zeros(1, round(par.grid_ioi*par.fs) - round(length(ir)))]; 
t_ir = [0: length(ir)-1] / par.fs; 
ax = pnl(size(par.all_pats, 1)+1, 5).select(); 
plot(ax, t_ir, ir, 'linew', 2, 'color', 'k'); 
ax.XLim = [0, 2.4]; 
ax.Visible = 'off';




ir = [ir, zeros(1, par.fs * 10)]; 

N = length(ir); 

mX = abs(fft(ir))  / N; 
freq = [0 : N-1] / N * par.fs; 

ax = pnl(size(par.all_pats, 1)+1, 6).select(); 
plot_fft(freq, mX, 'maxfreqlim', 15, 'linew', 0.1, 'ax', ax);
ax.YAxis.Visible = 'off';


%%

fname = '/datadisk/projects_backed_up/autocorrelation/figures/general/explain_5Hz/explain_5Hz'; 

print(fname, '-dsvg', '-painters', f);  
print(fname, '-dpng', '-painters', f);  

%%











