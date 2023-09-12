function main_coch(par)

tracks = {'H_unsyncopated', 'H_syncopated'}; 

responses = {'hilbert', 'slaney'}; 

n_tracks = length(tracks);

%%

    
for i_track=1:n_tracks

    track = tracks{i_track}; 

    col_names = {
        'response', 'z_meter_fft', 'z_meter_acf' ...
        };

    tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

    data_to_plot = []; 

    for i_resp=1:length(responses)

        
        if strcmp(responses{i_resp}, 'hilbert')
            [s, fs] = audioread(fullfile(...
                par.coch_data_path, ...
                sprintf(...
                '%s_200ms_calib75dB_max_short-term_ff_freq1237Hz_onset10_offset50.wav', track...
                ))); 
            s = s ./ max(abs(s)); 
            x = abs(hilbert(s')); 
            x = decimate(x, 10); 
            fs = fs / 10; 
        elseif strcmp(responses{i_resp}, 'slaney')
            coch = load(fullfile(...
                        par.coch_data_path, ...
                        'Slaney_128coch_meddis_timeDomain_meanAcrossCF.mat'...
                        )); 
            x = coch.hc(strcmp(coch.cond_names, track), :); 
            idx = round(0.010 * coch.fs); 
            x(1:idx) = x(idx); 
            x = decimate(x, 10); 
            fs = coch.fs / 10; 
        end

        t = [0:length(x)-1]/fs; 

        % withuout aperiodic subtraction    
        [acf, lags, ~, mX, freq] = get_acf(x, fs);    

        % get features 
        feat_acf = get_acf_features(acf, lags, ...
                                 par.lags_meter_rel, par.lags_meter_unrel);         

        feat_fft = get_fft_features(mX, freq, ...
                    par.freq_meter_rel, par.freq_meter_unrel); 
                

        new_row = [responses(i_resp), {feat_fft.z_meter_rel}, {feat_acf.z_meter_rel}]; 
        tbl = [tbl; new_row]; 

        data_to_plot(i_resp).response = responses{i_resp}; 
        data_to_plot(i_resp).x = x; 
        data_to_plot(i_resp).t = t; 
        data_to_plot(i_resp).mX = mX; 
        data_to_plot(i_resp).freq = freq; 
        data_to_plot(i_resp).acf = acf; 
        data_to_plot(i_resp).lags = lags; 


    end

    % save DATA 
    fname = sprintf('coch_track-%s', strrep(track, '_', ''));
    save(fullfile(par.data_path, [fname, '.mat']), 'data_to_plot', 'par'); 
    writetable(tbl, fullfile(par.data_path, [fname, '.csv'])); 



end


