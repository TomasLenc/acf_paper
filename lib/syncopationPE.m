function C = syncopationPE(in, grouping, varargin)
% 
% 
% 
% C = omega*O+U
% 
% counterevidence (C) against the induction of a meter is determined as 
% the weighted (omega) sum of omitted (O) and unaccented (U) downbeats
% 
% omega = 4:1 (weighting of silent vs. unaccented downbeats)
%

W = 4; 

% if any(strcmpi(varargin, 'perbar'))
%  
%     nbars = size(in,2)/grouping; 
%     C = zeros(1,size(in,1)); 
%     
%     for file=1:size(in,1)
%         pattern = in(file,:); 
%         accents = zeros(size(pattern)); 
%         for i=1:length(pattern)
%             if pattern(i)
%                 if i+1>length(pattern)
%                     accents(i)=1; 
%                     continue
%                 end
%                 if pattern(i+1)==0
%                     accents(i)=1; 
%                     continue
%                 end
%                 if (i-1<1 & pattern(i+1)==0) | (i-1<1 & pattern(i+1)==1 & pattern(i+2)==1)
%                     accents(i)=1; 
%                     continue
%                 end
%                 if (i-1<1 & pattern(i+1)==1 & pattern(i+2)==0)
%                     continue
%                 end
%                 if (pattern(i-1)==0 & pattern(i+1)==1 & i+2>length(pattern))
%                     continue
%                 end
%                 if (pattern(i-1)==0 & pattern(i+1)==1 & pattern(i+2)==1)
%                     accents(i)=1; 
%                     continue
%                 end
%             end
%         end
% 
%         coinc_accent = sum(accents(1:grouping:end)==1); 
%         coinc_unaccent = sum(accents(1:grouping:end)==0); 
%         coinc_silence = sum(pattern(1:grouping:end)==0); 
% 
%         C(file) = (W * coinc_silence) + (coinc_unaccent); 
%     end
%     
% 
%     
% else
%     
    
%     
%     C = zeros(1,size(in,1)); 
%     for file=1:size(in,1)
%         pattern = in(file,:); 
%         accents = zeros(size(pattern)); 
%         for i=1:length(pattern)
%             if pattern(i)
%                 if i+1>length(pattern)
%                     accents(i)=1; 
%                     continue
%                 end
%                 if pattern(i+1)==0
%                     accents(i)=1; 
%                     continue
%                 end
%                 if (i-1<1 & pattern(i+1)==0) | (i-1<1 & pattern(i+1)==1 & pattern(i+2)==1)
%                     accents(i)=1; 
%                     continue
%                 end
%                 if (i-1<1 & pattern(i+1)==1 & pattern(i+2)==0)
%                     continue
%                 end
%                 if (pattern(i-1)==0 & pattern(i+1)==1 & i+2>length(pattern))
%                     continue
%                 end
%                 if (pattern(i-1)==0 & pattern(i+1)==1 & pattern(i+2)==1)
%                     accents(i)=1; 
%                     continue
%                 end
%             end
%         end
% 
%         coinc_accent = sum(accents(1:grouping:end)==1); 
%         coinc_unaccent = sum(accents(1:grouping:end)==0); 
%         coinc_silence = sum(pattern(1:grouping:end)==0); 
% 
%         C(file) = (W * coinc_silence) + (coinc_unaccent); 
%     end
%     
%     
%     
    
    

% first pad the pattern from front and back as if it was repeated 
   
C = zeros(1,size(in,1)); 
for file=1:size(in,1)
    padded_pattern = repmat(in,1,3); 
    pattern = in(file,:); 
    accents = zeros(size(pattern)); 
    nevents = length(pattern); 
    for i=[1:length(pattern)]+nevents
        if padded_pattern(i)
            if padded_pattern(i+1)==0
                accents(i-nevents)=1; 
                continue
            end
            if (padded_pattern(i-1)==0 & padded_pattern(i+1)==1 & padded_pattern(i+2)==1)
                accents(i-nevents)=1; 
                continue
            end
        end
    end

    coinc_accent = sum(accents(1:grouping:end)==1); 
    coinc_unaccent = sum(pattern(1:grouping:end)==1 & accents(1:grouping:end)==0); 
    coinc_silence = sum(pattern(1:grouping:end)==0); 

    C(file) = (W * coinc_silence) + (coinc_unaccent); 
end



