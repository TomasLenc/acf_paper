function res = syncopationLHL(pattern, meter, events_in_cycle, varargin) 
% "If N is a note that precedes a rest, R, and R has a metric weight greater than or equal to N, 
% then the pair (N, R) is said to constitute a monophonic syncopation."
% 
% If N < R
% Syncopation = R - N

nbars = size(pattern, 2) / events_in_cycle; 

if strcmpi(meter,'[2,2]')
    salience = [0,-2,-1,-2]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
elseif strcmpi(meter,'[2,2,3]')
    salience = [0,-3,-2,-3, -1,-3,-2,-3, -1,-3,-2,-3]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
elseif strcmpi(meter,'[3,2,2]')
    salience = [0,-3,-3,-2,-3,-3, -1,-3,-3,-2,-3,-3]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
elseif strcmpi(meter,'[2,3,2]')
    salience = [0,-3,-2,-3,-2,-3, -1,-3,-2,-3,-2,-3]; 
    weightgrid = repmat(salience, 1, ceil(size(pattern,2)/length(salience))); 
else
    error('meter not implemented')
end
    
    
if any(strcmpi(varargin, 'perbar'))
    res = zeros(size(pattern,1),nbars); 
    
    for i_event=1:size(pattern,1)
        synidx = 0; 
        bar=1; 
        for i_event=1:size(pattern,2)
            if pattern(file,i_event)
                c=1; 
                tmpx = []; 
                while 1
                    if (i_event+c)>size(pattern,2)
                        break
                    elseif pattern(file,i_event+c)==1
                        break
                    elseif pattern(file,i_event+c)==0
                        tmpx = [tmpx, weightgrid(i_event+c)]; 
                    end
                    c = c+1; 
                end
                if any(tmpx>weightgrid(i_event))
                    synidx = synidx + (max(tmpx)-weightgrid(i_event)); 
                end
            end
            if mod(i_event+1,events_in_cycle)==0
               res(file, bar) = synidx; 
               synidx=0; 
               bar = bar+1;  
            end
        end
    end
    
else
    
    res = zeros(1,size(pattern,1)); 

    for file=1:size(pattern,1)
        synidx = 0; 
        for i_event=1:size(pattern,2)
            if pattern(file,i_event)
                c=1; 
                tmpx = []; 
                while 1
                    if (i_event+c)>size(pattern,2)
                        break
                    elseif pattern(file,i_event+c)==1
                        break
                    elseif pattern(file,i_event+c)==0
                        tmpx = [tmpx, weightgrid(i_event+c)]; 
                    end
                    c = c+1; 
                end
                if any(tmpx>weightgrid(i_event))
                    synidx = synidx + (max(tmpx)-weightgrid(i_event)); 
                end
            end
        end

        res(file) = synidx; 
    end

    
end



