function vts = get_vowel_t(pa, vidx, varargin)
nv = length(vidx);
vts = nan(1, nv);

if ~isempty(fsic(varargin, 'peakRMS'))
    vtMode = 'peakRMS';
    
    t_rms = varargin{fsic(varargin, 'peakRMS') + 1};
    t0 = varargin{fsic(varargin, 'peakRMS') + 2};
else
    vtMode = 'tMid';
end

if isequal(vtMode, 'peakRMS')
%     TODO;

%     for i1 = 1 : nv
%         tv0 = pa.tbeg(vidx(i1));
%         tv1 = pa.tend(vidx(i1));
%         
%         idx0 = round(tv0 / t0);
%         idx1 = round(tv1 / t0);
%         v_rms = t_rms(idx0 : idx1);
%     end
else
    for i1 = 1 : nv
        tv0 = pa.tbeg(vidx(i1));
        tv1 = pa.tend(vidx(i1));
        
        vts(i1) = (tv0 + tv1) / 2;
    end
end
return