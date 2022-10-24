function feat = get_acf_features(acf, lags_time, ...
                                 lags_meter_rel, lags_meter_unrel, ...
                                 lags_meter_unrel_left, lags_meter_unrel_right)
% Extract features from acf. Time lag must be last dimension. 

% get indices for lags of interest
lags_meter_rel_idx = dsearchn(lags_time', lags_meter_rel')'; 
lags_meter_unrel_idx = dsearchn(lags_time', lags_meter_unrel')'; 
lags_meter_unrel_left_idx = dsearchn(lags_time', lags_meter_unrel_left')'; 
lags_meter_unrel_right_idx = dsearchn(lags_time', lags_meter_unrel_right')'; 


% calculate mean acf value for lags of interest
subs_cmd = []; 
subs_cmd.subs = repmat({':'}, 1, ndims(acf)); 
subs_cmd.type = '()'; 

subs_cmd.subs{end} = lags_meter_rel_idx;
acf_mean_meter_rel = mean(subsref(acf, subs_cmd), ndims(acf)); 

subs_cmd.subs{end} = lags_meter_unrel_idx;
acf_mean_meter_unrel = mean(subsref(acf, subs_cmd), ndims(acf)); 

subs_cmd.subs{end} = lags_meter_unrel_left_idx;
acf_mean_meter_unrel_left = mean(subsref(acf, subs_cmd), ndims(acf)); 

subs_cmd.subs{end} = lags_meter_unrel_right_idx;
acf_mean_meter_unrel_right = mean(subsref(acf, subs_cmd), ndims(acf)); 

% z-score
subs_cmd.subs{end} = [lags_meter_rel_idx, lags_meter_unrel_idx];
z = zscore(subsref(acf, subs_cmd), [], ndims(acf)); 

subs_cmd.subs{end} = [1 : length(acf_mean_meter_rel)];
z_meter_rel = mean(subsref(z, subs_cmd), ndims(acf)); 

% acf ratio 
ratio_meter_rel = acf_mean_meter_rel ./ acf_mean_meter_unrel;
ratio_meter_rel_left = acf_mean_meter_rel ./ acf_mean_meter_unrel_left;
ratio_meter_rel_right = acf_mean_meter_rel ./ acf_mean_meter_unrel_right;

% acf contrast 
contrast_meter_rel = (acf_mean_meter_rel-acf_mean_meter_unrel) ./ ...
                     (acf_mean_meter_rel+acf_mean_meter_unrel); 

feat = []; 
feat.z_meter_rel = z_meter_rel; 
feat.ratio_meter_rel = ratio_meter_rel; 
feat.ratio_meter_rel_left = ratio_meter_rel_left; 
feat.ratio_meter_rel_right = ratio_meter_rel_right; 
feat.contrast_meter_rel = contrast_meter_rel; 



