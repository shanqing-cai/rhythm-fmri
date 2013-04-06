function [p1, p2, p3] = get_triphone_parts(mt)
idx_min = strfind(mt, '-');
idx_plu = strfind(mt, '+');

if length(idx_min) == 1 && length(idx_plu) == 1
    p1 = mt(1 : idx_min - 1);
    p2 = mt(idx_min + 1 : idx_plu - 1);
    p3 = mt(idx_plu + 1 : end);
elseif length(idx_min) == 1
    p1 = mt(1 : idx_min - 1);
    p2 = mt(idx_min + 1 : end);
    p3 = '';
elseif length(idx_plu) == 1
    p1 = '';
    p2 = mt(1 : idx_plu - 1);
    p3 = mt(idx_plu + 1 : end);
else
    p1 = '';
    p2 = mt;
    p3 = '';
end
return