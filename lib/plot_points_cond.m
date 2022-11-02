function plot_points_cond(ax, x1, x2, col1, col2)
% plot two sets of scalar values as points to compare the two groups

x = {x1, x2}; 
cols = {col1, col2}; 

hold(ax, 'on')

for i_cond=1:length(x)
%     plot(ax, i_cond, x{i_cond}, 'o', 'color', cols{i_cond})
%     mu = mean(x{i_cond});
%     ci = std(x{i_cond}) / sqrt(length(x{i_cond})) * norminv(0.95); 
%     plot(ax, i_cond, mu, 'o', 'color', 'k', 'MarkerFaceColor', 'k')
%     plot(ax, [i_cond, i_cond], [mu-ci, mu+ci], 'color', 'k', 'linew', 1)
     plot_points(ax, i_cond, x{i_cond}, 'col', cols{i_cond})
end

ax.XLim = [0.5, 2.5]; 

[~, p] = ttest(x1, x2); 

title(ax, sprintf('p=%.5f', p)); 