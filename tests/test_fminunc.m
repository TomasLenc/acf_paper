
theta = [5, 0.5, 2]; 

x = [0.1 : 0.1 : 100]; 
y = aperiodic_knee(theta, x); 

theta_opt = fit_aperiodic(x, y); 

theta
theta_opt

y_pred = aperiodic_knee(theta_opt, x); 

figure
plot(x, y, 'k', 'linew', 2)
hold on 
plot(x, y_pred, 'r--', 'linew', 2)