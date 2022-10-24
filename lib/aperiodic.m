function y = aperiodic(theta, x, knee)

if knee
    y = theta(1) - log(theta(2) + x .^ theta(3)); 
else
    y = theta(1) - log(x .^ theta(end)); 
end
