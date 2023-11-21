function s_old = appendAllFields(s_old, s_new)

names = fieldnames(s_old);

for i_name=1:length(names)
    
    shape = size(s_old.(names{i_name})); 
    s_old.(names{i_name}) = cat(1, s_old.(names{i_name}), s_new.(names{i_name})); 
    
end