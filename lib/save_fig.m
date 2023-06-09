function save_fig(f, fname)

% strip extension
if ismember(fname(end-3:end), {'.png', '.eps', '.svg', '.fig'})
    fname = fname(1:end-4);
end

% saveas(f, fullfile(par.fig_path, [fname, '.fig']));  
% print(fullfile(par.fig_path, [fname, '.png']), '-dpng', '-painters', f);  
print(fullfile([fname, '.svg']), '-dsvg', '-painters', f);  

