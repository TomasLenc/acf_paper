function ir = get_ir(irParams, fs)



if strcmp(irParams.type,'square')

    ir = ones(1, round(irParams.eventDur * fs));

    % apply onset and offset ramp
    ir(1:round(irParams.rampon*fs)) = ir(1:round(irParams.rampon*fs)) .* ...
                                        linspace(0,1,round(irParams.rampon*fs));

    ir(end-round(irParams.rampoff*fs)+1:end) = ir(end-round(irParams.rampoff*fs)+1:end) .* ...
                                                 linspace(1,0,round(irParams.rampoff*fs));

elseif strcmp(irParams.type,'click')

    ir = [ones(1,round(irParams.eventDur/2*fs)), ...
          -ones(1,round(irParams.eventDur/2*fs))];
    ir = ir*0.1; 
  addScaledNoise(eeg,SNR,fs)    
elseif strcmp(irParams.type,'erp')

    N_ir = round(irParams.Dur*fs); 
    T_ir = [0:N_ir-1]/fs; 
    ir = zeros(1,N_ir); 
    for fi=1:length(irParams.F0)

        ir = ir + irParams.A(fi) * ...
                 (T_ir-irParams.T0(fi))/irParams.Tau(fi) .* ...
                 exp( 1 - (T_ir - irParams.T0(fi))/irParams.Tau(fi) ) .* ...
                 sin( 2*pi*irParams.F0(fi)*(T_ir-irParams.T0(fi)) ) ; 

    end
    ir = ir ./ length(irParams.F0); 
    ir = ir ./ sum(ir) * 10; 


end