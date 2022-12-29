function par = get_par()

acf_tools_path = '/datadisk/projects_git_dl/acf_tools/src'; 
rnb_tools_path = '/datadisk/projects_git_dl/rnb_tools/src'; 

lw_path = '/datadisk/projects_backed_up/autocorrelation/simulations/lib_external/letswave6'; 
pica_path = '/datadisk/projects_backed_up/autocorrelation/simulations/lib_external/piCA'; 

%% return structure 

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end


