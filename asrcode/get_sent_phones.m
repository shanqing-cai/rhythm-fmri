function phns = get_sent_phones(stc, stcs, sentPhns)
idx = fsic(stcs, stc);

if isempty(idx)
    phns = [];
else
    phns = sentPhns{idx};
end

return