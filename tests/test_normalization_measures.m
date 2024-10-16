
idx_rel = [1, 2, 3]; 
idx_unrel = [4:100]; 

vals = rand(1, 100); 

%% 

scale_factors = [1 : 100]; 

mu = []; 
sigma = []; 
rat = []; 
norm_sd = []; 

for i=1:length(scale_factors)

    x = vals * scale_factors(i); 
    
    mu(i) = mean(x); 
    sigma(i) = std(x); 
    
    norm_mean(i) = (mean(x(idx_rel)) - mean(x(idx_unrel))) / mean(x); 
    norm_sd(i) =  (mean(x(idx_rel)) - mean(x(idx_unrel))) / std(x); 
end

figure
plot(mu, sigma)

figure
plot(norm_mean)
figure
plot(norm_sd)

%%

