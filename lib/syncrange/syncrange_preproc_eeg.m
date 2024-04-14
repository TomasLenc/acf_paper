clear 

par = get_par();

data_path = '/DATA1/XPSyncRange/individual_data'; 

tbl_subj = readtable('xpsyncrange_preproc_dirnames.txt'); 

trial_duration = 40.8; 

rhythms = {
        '6' 
        '42'
        '37'
        '31'
        '17'
        '26'
        '41'
        '19'
    };

n_trials_all = [9:-1:1]; 

%%

n_sub = size(tbl_subj, 1);
n_rhythms = length(rhythms); 

for i_n_trials=1:length(n_trials_all)
    
    n_trials = n_trials_all(i_n_trials); 

    for i_rhythm=1:n_rhythms
    
        rhythm_id = rhythms{i_rhythm};
        
        data_all = {1, n_sub}; 
        
        for i_sub=1:n_sub
            
            sub_id = tbl_subj{i_sub, 1}{:};
            sub_preproc_dir = tbl_subj{i_sub, 2}{:};
            
            fprintf('rhythm-%s_nTrials-%d_sub-%s\n', rhythm_id, n_trials, sub_id);
    
            % load log
            d_log = dir(fullfile(data_path, sub_id, 'log*.mat'));
            log = load(fullfile(d_log.folder, d_log.name));
            
            % load preprocessed data
            fpath = fullfile(data_path, sub_id, 'lw', sub_preproc_dir); 
            fname = sprintf('rhythm%s *.lw6', rhythm_id); 
            d = dir(fullfile(fpath, fname));
            
            [header, data] = CLW_load(fullfile(d.folder, d.name));
            
            [header, data] = RLW_arrange_epochs(header, data, ...
                                    [1:min(n_trials, header.datasize(1))]); 

            [header, data] = RLW_average_epochs(header, data); 
            
            [header, data] = RLW_rereference(header, data, ...
                                    'apply_list', {header.chanlocs.labels}, ...
                                    'reference_list',{'mast1', 'mast2'}); 
            
            [header, data] = RLW_arrange_channels(header, data,  ...
                                    {'F1','Fz','F2','FC1','FCz','FC2','C1','Cz','C2',}); 
            
            [header, data] = RLW_butterworth_filter(header, data, ...
                                                    'filter_type', 'lowpass', ...
                                                    'high_cutoff', 20, ...
                                                    'filter_order', 2); 
                                                
            [header, data] = segment_safe(header, data, {'1'}, ...
                                          'x_start', 0, 'x_duration', trial_duration, ...
                                          'ignore_out_of_range', true); 
                                      
            [header, data] = RLW_pool_channels(header, data, {header.chanlocs.labels},...
                                                'keep_original_channels', false);
            
            data_all{i_sub} = squeeze(data)'; 
            
        end
        
        % save
        fs = 1/header.xstep; 
        
        data = cat(1, data_all{:});
        
        save_fpath = fullfile(par.eeg_path, 'syncrange');
        save_fname = sprintf('exp-syncrange_rhythm-%s_nTrials-%d_eeg.mat', ...
                             rhythm_id, n_trials);
        
        if ~isdir(save_fpath)
            mkdir(save_fpath)
        end
    
        save(fullfile(save_fpath, save_fname), 'data', 'fs');
    
    end

end



%% erp cycle plot
% 
% x = epoch_chunks(data, fs, 2.4);
% x = squeeze(mean(x, 1)); 
% t = [0 : length(x)-1]/fs;
% 
% idx = find(strcmp(log.res.blockOrder, rhythm_id));
% s_fs = log.res.data.block(idx).trial(1).stimulus.fs; 
% s = log.res.data.block(idx).trial(1).stimulus.s; 
% s = epoch_chunks(s, s_fs, 2.4);
% s = mean(s, 1);  
% s = (s - min(s)) ./ (max(s)-min(s)) ;
% s = s * (max(x) - min(x)) + min(x); 
% t_s = [0 : length(s)-1]/s_fs; 
% 
% figure('color', 'white', 'Position', [675 735 363 227])
% plot(t_s, s, 'linew', 2, 'color', [.7, .7, .7]); 
% hold on
% plot(t, x, 'linew', 2); 
% box off
% 










