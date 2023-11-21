% FFT of ACF as IR changes

clear 

par = get_par(); 

par.duty_cycles = linspace(0.050, 0.180, 6); 

f = figure('color', 'white', 'Position', [651 366 532 574]); 
pnl = panel(f); 
pnl.pack('h', 4);
for i=1:4
    pnl(i).pack('v', 6); 
end

for i_cond=1:6
    
    ir = get_square_kernel(par.fs, ...
        'duration', par.duty_cycles(i_cond), ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 

    % make whole signal 
    [x, t] = get_s(...
                    par.pat, ...
                    par.grid_ioi, ...
                    par.fs, ...
                    'n_cycles', par.n_cycles, ...
                    'ir', ir ...
                    );

    % get acf withuout aperiodic subtraction    
    [acf, lags, ~, mX, freq] = get_acf(x, par.fs);    

    mX_acf = abs(fft(acf)); 
    freq_acf = [0 : length(mX_acf)-1] /  length(mX_acf) * par.fs; 

    feat_fft = get_fft_features(mX, freq, ...
                            par.freq_meter_rel, par.freq_meter_unrel); 
    feat_fft_acf = get_fft_features(mX_acf, freq_acf, ...
                            par.freq_meter_rel, par.freq_meter_unrel); 

    ax = pnl(1, i_cond).select(); 
    plot(t, x); 
    xlim([0, 2.4])
    
    ax = pnl(2, i_cond).select(); 
    plot(freq, mX)
    xlim([0, 4.9])
    title(sprintf('z = %.3f', feat_fft.z_meter_rel)); 

    ax = pnl(3, i_cond).select(); 
    plot(lags, acf)
    xlim([0, 4.8])
    
    ax = pnl(4, i_cond).select(); 
    plot(freq_acf, mX_acf)
    xlim([0, 4.9])
    title(sprintf('z = %.3f', feat_fft_acf.z_meter_rel)); 
    
    if i_cond<6
        for i=1:4
            ax = pnl(i, i_cond).select(); 
            ax.XTick = []; 
        end
    end
    
end

pnl.de.margin = 1;
pnl.de.marginleft = 15;
pnl.de.margintop = 10;
pnl.margintop = 15; 

i=1;
tit = pnl(i).title('time'); 
tit.Position(2) = tit.Position(2) * 1.07; 

i=2;
tit = pnl(i).title('FFT'); 
tit.Position(2) = tit.Position(2) * 1.07; 

i=3;
tit = pnl(i).title('ACF'); 
tit.Position(2) = tit.Position(2) * 1.07; 

i=4; 
tit = pnl(i).title('FFT of ACF'); 
tit.Position(2) = tit.Position(2) * 1.07; 


saveas(f, '/home/tomo/Downloads/fft_of_acf_duty.fig'); 
saveas(f, '/home/tomo/Downloads/fft_of_acf_duty.png'); 