function all_good_pats = find_all_patterns(n_events, n_sounds, max_group_size)
% finds all possible patterns that have n_sounds in n_events, rejecting all
% mirrored versions, circshifted versions, and allowing up to max_group_size
% successive sounds. 

f_gen_all_pats = @(n, k) dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0'; 
all_possible_pats = f_gen_all_pats(n_events, n_sounds); 
all_possible_pats = flip(all_possible_pats, 1); 

all_good_pats = nan(size(all_possible_pats)); 

c = 1; 
for i_pat=1:size(all_possible_pats, 1)
       
    pat = all_possible_pats(i_pat, :); 

    flag = 1; 
    
    for i_shift=0:n_events
        [B, N] = RunLength(circshift(pat, i_shift));
        if any(N(B == 1) > max_group_size) 
            flag = 0; 
        end
    end
    if flag == 0
        continue
    end
    
    for i_shift=0:n_events
        if any(all((all_good_pats - circshift(pat, i_shift)) == 0, 2)) || ...
           any(all((all_good_pats - circshift(flip(pat), i_shift)) == 0, 2))
            flag = 0; 
        end
    end    
    if flag == 0
        continue
    end
   
    all_good_pats(c, :) = pat; 
    
    c = c+1; 
end

all_good_pats(c:end, :) = []; 
