function lags = get_lag_harmonics(lag, max_lag, varargin)

parser = inputParser; 

addParameter(parser, 'lag_harm_to_exclude', []); 

parse(parser, varargin{:}); 

lag_harm_to_exclude = parser.Results.lag_harm_to_exclude; 


lags = [lag : lag : max_lag]; 

mask_to_rm = zeros(1, length(lags), 'like', false); 

for i_excl=1:length(lag_harm_to_exclude)
    
    lag_to_exclude_harm = [...
        lag_harm_to_exclude(i_excl) : lag_harm_to_exclude(i_excl) : max_lag...
        ]; 
        
    for li=1:length(lags)
       
        if any(abs(lags(li) - lag_to_exclude_harm) < 1e-8)
            mask_to_rm(li) = true; 
        end
        
    end
    
end

lags = lags(~mask_to_rm); 