function mvl = get_mean_vwl_level(data, asrPAlign, vidx, micRMS_100dBA)
vwlInts = {};

if numel(asrPAlign.phones) <= 2
    mvl = NaN;
    return
end

% vidx = get_vowel_indices(data.params.name);

frameLen = data.params.frameLen / data.params.sr;
meanRMS_vwls = nan(1, length(vidx));
for i1 = 1 : numel(vidx)   
    vwlInts{end + 1} = [round(asrPAlign.tbeg(vidx(i1) + 1) / frameLen), ...
                        round(asrPAlign.tend(vidx(i1) + 1) / frameLen)];

    if vwlInts{end}(2) > vwlInts{end}(1)
        t_rms = data.rms(vwlInts{end}(1) : vwlInts{end}(2), 1);
        meanRMS_vwls(i1) = rms(t_rms);    
    end
end

mvl = rms(meanRMS_vwls(~isnan(meanRMS_vwls)));
mvl = 100 + 20 * log10(mvl / micRMS_100dBA);

%% 
return