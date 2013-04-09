function figIdDat=makeFigDataMon()
fid=figure('Position', [300, 100, 1200, 640], ...
           'Name', 'Data monitor', 'NumberTitle', 'off');
axes1=subplot('Position', [0.05 ,0.575, 0.25, 0.4]);   % Input waveform
axes2=subplot('Position', [0.05, 0.1, 0.25, 0.4]);   % Output waveform

axes4=subplot('Position', [0.35, 0.23, 0.30, 0.14]);
axes5=subplot('Position', [0.35, 0.05, 0.64, 0.14]);
axes7=subplot('Position', [0.69, 0.23, 0.30, 0.14]);

axes3 = subplot('Position', [0.35, 0.650, 0.64, 0.100]);
axes6 = subplot('Position', [0.35, 0.750, 0.64, 0.200]);
axes8 = subplot('Position', [0.35, 0.435, 0.64, 0.200]);

% axes5=subplot('Position',[0.7,0.575,0.275,0.4]);
% axes6=subplot('Position',[0.7,0.1,0.275,0.4]);
figIdDat = [fid, axes1, axes2, axes3, axes4, axes5, axes6, axes7, axes8];
return