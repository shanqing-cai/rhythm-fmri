function phns = get_dict_pronun(word, dictFN)
word = upper(word);

if ~iscell(dictFN)
    dtxt = textread(dictFN, '%s', 'delimiter', '\n');
else
    dtxt = dictFN;
end

phns = {};
for i1 = 1 : length(dtxt)
    t_line = dtxt{i1};
    if length(t_line) > length(word)
        if isequal(t_line(1 : length(word)), word)
            t_items = splitstring(t_line);
            for j1 = 2 : length(t_items)
                phns{end + 1} = strrep(strrep(strrep(t_items{j1}, '0', ''), '1', ''), '2', '');
            end
            break;
        end
    end
end
return