function [s_trial, t] = get_s(pattern, n_cycles, grid_ioi, jitter, ir, fs, ...
                              varargin)

% magnitude of periodic emphasis as a proportion of max amplitude                          
emph_magn = 0; 
if any(strcmp(varargin, 'emph_magn'))
    emph_magn = varargin{find(strcmp(varargin, 'emph_magn')) + 1}; 
end
% period of periodic emphasis (in units of grid points)                  
emph_period = 4; 
if any(strcmp(varargin, 'emph_period'))
    emph_period = varargin{find(strcmp(varargin, 'emph_period')) + 1}; 
end
% phase of periodic emphasis (in units of grid points)                          
emph_phase = 0; 
if any(strcmp(varargin, 'emph_phase'))
    emph_phase = varargin{find(strcmp(varargin, 'emph_phase')) + 1}; 
end

                         
% calculate duration of one long trial in seconds
trial_dur = n_cycles * length(pattern) * grid_ioi;

% make time vector for one trial
t = [0 : 1/fs : trial_dur-1/fs];

% allocate envelope vector for the whole trial (as zeros)
s_trial = zeros(size(t));

% make ncycles copies of the rhythmic pattern
pattern_whole_trial = repmat(pattern, 1, n_cycles);

emph_magn_abs = max(abs(pattern_whole_trial)) * emph_magn; 

% go over each event in the trial
for i=1:length(pattern_whole_trial)
    
    % find the time of event onset
    event_time = (i-1)*grid_ioi;
    % apply jitter
    event_time = event_time + jitter*randn(); 
    event_time = max(event_time, 0); 
    event_time = min(event_time,trial_dur); 
    % convert to index
    event_idx = round(event_time*fs);
    
    amp = pattern_whole_trial(i); 
    % find whether there's emphasis on this grid point
    if mod((i-1-emph_phase), emph_period) == 0
        amp = amp + emph_magn_abs; 
    end

    s_trial(event_idx+1) = amp;
end

% convolve with impulse response
s_trial = conv(ir, s_trial, 'full'); 
s_trial = s_trial(1:length(t)); 