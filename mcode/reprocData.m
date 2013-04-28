function dataOut = reprocData(data, varargin)
p=data.params;

if isempty(fsic(varargin, 'sig'))
    sigIn = data.signalIn;
else
    sigIn = varargin{fsic(varargin, 'sig') + 1};
end

if ~isempty(varargin)
    for i1=1:2:length(varargin)
        if ~isequal(varargin{i1}, 'sig')
            p.(varargin{i1}) = varargin{i1 + 1};
        end
    end
end
    
MexIO('reset');
MexIO('init', p);

sigIn=resample(sigIn, data.params.sr * data.params.downfact, data.params.sr);
sigInCell=makecell(sigIn, data.params.frameLen * data.params.downfact);
for n = 1 : length(sigInCell)
%     tic;
    TransShiftMex(5,sigInCell{n});
end


dataOut=MexIO('getData');



return