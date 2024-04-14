clear 

par = get_par();

data_path = '/DATA1/XPSyncRange/tapping/individual_data'; 

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

%%

n_sub = size(tbl_subj, 1);
n_rhythms = length(rhythms); 

for i_rhythm=1:n_rhythms

    rhythm_id = rhythms{i_rhythm};
    
    data_all = {1, n_sub}; 
    
    for i_sub=1:n_sub
        
        sub_id = tbl_subj{i_sub, 1}{:};

        fprintf('rhythm-%s_sub-%s\n', rhythm_id, sub_id);

        sub_raw_dir = fullfile(data_path, sub_id);
        
        % load log
        d_log = dir(fullfile(sub_raw_dir, 'log*.mat'));
        log = load(fullfile(d_log.folder, d_log.name));
        
        % load raw data
        fname = sprintf('*_rhythmID%s.wav', rhythm_id); 
        d = dir(fullfile(sub_raw_dir, fname));
        
        [x, fs] = audioread(fullfile(d.folder, d.name));

        N = round(trial_duration * fs); 
        x = x(1:N); 

        % downsample
        ds_ratio = 80;
        fs_ds = fs / ds_ratio;
        x_ds = decimate(x, ds_ratio);

        data_all{i_sub} = x_ds'; 
        
    end
    
    % save    
    data = cat(1, data_all{:});
    
    save_fpath = fullfile(par.eeg_path, 'syncrange');
    save_fname = sprintf('exp-syncrange_rhythm-%s_tapping.mat', ...
                         rhythm_id);
    
    if ~isdir(save_fpath)
        mkdir(save_fpath)
    end

    fs = fs_ds; 
    save(fullfile(save_fpath, save_fname), 'data', 'fs');

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










