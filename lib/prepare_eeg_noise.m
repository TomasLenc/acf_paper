function noise = prepare_eeg_noise(n_rep, trial_dur)

par = get_par(); 

% get a list of all available subjects 
d_subjects = dir(fullfile(par.resting_eeg_path, 'sub-*'));

% allocate 
noise = cell(1, n_rep);

trial_dur = par.n_cycles * length(par.pat) * par.grid_ioi; 

parfor i_rep=1:n_rep

    % randomly pick a subject 
    i_sub = randsample(length(d_subjects), 1);

    d_eeg = dir(fullfile(d_subjects(i_sub).folder, d_subjects(i_sub).name,  ...
        'ses-session1', 'eeg', ...
        'sub-*_ses-session1_task-eyesopen_eeg.vhdr'));

    % load EEG data 
    [header, data] = RLW_import_VHDR(...
        fullfile(d_eeg.folder, d_eeg.name)); 

    % take average reference 
    [header, data] = RLW_rereference(header, data,...
        'apply_list', {header.chanlocs.labels}, ...
        'reference_list', {header.chanlocs.labels});

    % pick a random channel                            
    chan_idx = randsample(header.datasize(2), 1);

    [header, data] = RLW_arrange_channels(header, data, ...
                            {header.chanlocs(chan_idx).labels});

    % resample to match fs        
    [P, Q] = rat(par.fs / (1/header.xstep));
    X = squeeze(data)';
    X = resample(X, P, Q);
    data_resampled = [];
    data_resampled(1, 1, 1, 1, 1, :) = X';
    header_resampled = header;
    header_resampled.xstep = 1/par.fs;
    header_resampled.datasize = size(data_resampled);        

    % pick a starting point at random
    eeg_trial_dur = header.xstep * header.datasize(end) - header.xstart;       
    x_start = rand * (eeg_trial_dur - trial_dur);

    % segment the duration we need
    [header_ep, data_ep] = RLW_crop(header_resampled, data_resampled, ...
                    'x_crop', true, ...
                    'x_start', x_start, ...
                    'x_size', round(trial_dur / header_resampled.xstep) ...
                    );

    % save the data
    noise{i_rep} = squeeze(data_ep);                     

    fprintf('%d/%d added %s chan-%d from %.1fs \n', ...
        i_rep, n_rep, d_subjects(i_sub).name, chan_idx, x_start ...
        );

end

noise = cat(2, noise{:})'; 

% mean-center the data         
noise = noise - mean(noise, 2);       
