function format_pcf(base_pcf, out_pcf, ...
                    timeWarp_initStat, ...
                    timeWarp_tBegin, timeWarp_rate1, timeWarp_dur1, ...
                    timeWarp_durHold, timeWarp_rate2, ...
                    F1ShiftRatio, FmtShiftStat0, FmtShiftStat1)
%% Input sanity check
if FmtShiftStat1 < FmtShiftStat0
    error('FmtShiftStat1 must be >= FmtShiftStat0');
end

%%
base_txt = textread(base_pcf, '%s', 'delimiter', '\n');

f = fopen(out_pcf, 'wt');
inSect2 = 0;
for i1 = 1 : numel(base_txt)
    if i1 ~= 3
        if ~inSect2
            bDirectCopy = 1;
        else
            if length(strfind(base_txt{i1}, ', ')) == 4
                t_items = strsplit(base_txt{i1}, ', ');
                t_stat = str2double(t_items{1});
                if t_stat >= FmtShiftStat0 && t_stat <= FmtShiftStat1
                    bDirectCopy = 0;
                    newLine = sprintf('%d, 0, 0, %.3f, 0', t_stat, F1ShiftRatio);
                else
                    bDirectCopy = 1;
                end
            else
                bDirectCopy = 1;
            end
        end
        
        if bDirectCopy
            fprintf(f, '%s\n', base_txt{i1});
        else
            fprintf(f, '%s\n', newLine);
        end
    else
        fprintf(f, '%d, %.3f, %.3f, %.3f, %.3f, %.3f\n', ...
                    timeWarp_initStat, timeWarp_tBegin, timeWarp_rate1, ...
                    timeWarp_dur1, timeWarp_durHold, timeWarp_rate2);
    end
    
    if ~isempty(strfind(base_txt{i1}, 'Section 2:'))
        inSect2 = 1;
    end
end

fclose(f);

return