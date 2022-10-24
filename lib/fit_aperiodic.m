function theta_opt = fit_aperiodic(x, y, knee)

options = optimoptions('fminunc', ...
                       'MaxFunctionEvaluations', 1000, ...
                       'Display', 'off'); 

if nargin == 3 && knee
    theta_init = [1, 1, 1]; 
else
    knee = false; 
    theta_init = [1, 1]; 
end

theta_opt = fminunc(@(theta) sum((y - aperiodic(theta, x, knee)).^2), ...
                    theta_init, options); 
