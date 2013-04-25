function [ostMat, timeWarpCfg, rmsThresh, maxInterval_5_10] = ...
            get_ost_pcf_params_rhy(varargin)
% Syntax:
%       get_ost_params_rhy(); %  For testing only 
%       get_ost_params_rhy(inputDir, sent, bPlot, opts);
%           Opts: -t: run test perturbations
%                 warpOnsetTime:    default: 0.0; unit: s
%                 accelWarpRate:    default: 0.25
%                 decelWarpRate:    default: 2.00
%                 F1ShiftRatio:     default: 0.25
%                 FmtShiftStat0:    default: 5
%                 FmtShiftStat1:    default: 9
%                 basePCF:

%% Input
% inputDir = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_1\pre';
if nargin == 0
    inputDir = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_3\pre';

    sent = 'The steady bat gave birth to pups';

    bPlot = 0;
    
    bRunTest = 1;       
elseif nargin >= 3
    inputDir = varargin{1};
    sent = varargin{2};
    bPlot = varargin{3};
    
    bRunTest = ~isempty(fsic(varargin, '-t'));
end



%% CONFIG
ASR_CODE_PATH = 'e:\speechres\rhythm-fmri\asrcode\';
typeColors.R = [0, 0, 1];
typeColors.N = [0, 0.5, 0];

DEFAULT_RMS_SLOPE_WIN = 0.030000;

basePCF = 'E:\speechres\pip\pipcfg\rhy_pert_test_1.pcf';

warpOnsetTime = 0.0;
decelWarpRate = 0.25;
accelWarpRate = 2.0;

F1ShiftRatio = 0.25;
FmtShiftStat0 = 5;
FmtShiftStat1 = 9;

%% Additional input options
if ~isempty(fsic(varargin, 'warpOnsetTime'))
    warpOnsetTime = varargin{fsic(varargin, 'warpOnsetTime') + 1};
end

if ~isempty(fsic(varargin, 'decelWarpRate'))
    decelWarpRate = varargin{fsic(varargin, 'decelWarpRate') + 1};
end

if ~isempty(fsic(varargin, 'accelWarpRate'))
    accelWarpRate = varargin{fsic(varargin, 'accelWarpRate') + 1};
end

if ~isempty(fsic(varargin, 'F1ShiftRatio'))
    F1ShiftRatio = varargin{fsic(varargin, 'F1ShiftRatio') + 1};
end

if ~isempty(fsic(varargin, 'FmtShiftStat0'))
    FmtShiftStat0 = varargin{fsic(varargin, 'FmtShiftStat0') + 1};
end

if ~isempty(fsic(varargin, 'FmtShiftStat1'))
    FmtShiftStat1 = varargin{fsic(varargin, 'FmtShiftStat1') + 1};
end

%% Phones look up
if isequal(sent, 'The steady bat gave birth to pups');
    sentPhns = {'dh', 'ah', 's',  't',  'eh', ...
                'd',  'iy', 'b',  'ae', 't', ...
                'g',  'ey', 'v',  'b',  'er', ...
                'th', 't',  'uw', 'p',  'ah', ...
                'p',  's'};
else
    error('Unsupported sentence: %s', sent)
end
nPhones = length(sentPhns) + 2; % Includes the slis

%% Initialize variables for data holding
if isequal(sent, 'The steady bat gave birth to pups')
    protoStruct = struct('N', [], 'R', []);
    
    min_rms_dh_ah = protoStruct;
    max_rms_dh_ah = protoStruct;
    
    min_rms_ratio_ah = protoStruct;
    max_rms_ratio_s_t = protoStruct;
    
    min_rms_s_t_eh = protoStruct;
    max_rms_eh = protoStruct;
    
    max_d_iy_neg_str_len = protoStruct;
    max_eh_d_neg_str_len = protoStruct;
    
    max_d_iy_neg_str_span = protoStruct;
    max_eh_d_neg_str_span = protoStruct;
    
    max_d_iy_pos_str_len = protoStruct;
    
    min_rms_s2d = protoStruct;
    max_rms_s2d = protoStruct;
    
    stat_0_prm1 = protoStruct;
    stat_0_prm2 = protoStruct;
    
    stat_2_prm1 = protoStruct;
    stat_2_prm2 = protoStruct;
    
    steady_s_dur = protoStruct;
    
%     len_s_t_eh_d = protoStruct;
    len_eh_d = protoStruct;
else
    error('The sentence "%s" is currently not supported', sent);
end


%% Load data
if bPlot
    figure;
    hold on;
end

if ~iscell(inputDir)
    inputDirs = {inputDir};
else
    inputDirs = inputDir;
end

for i0 = 1 : numel(inputDirs)
    check_dir(inputDirs{i0});
end

for i0 = 1 : numel(inputDirs)
    inputDir = inputDirs{i0};
    
    drep = dir(fullfile(inputDir, 'rep*'));
    if isempty(drep)
        fprintf(2, 'WARNING: cannot find matching reps in directory: %s\n', inputDir);
        return;
    end
    for i1 = 1 : numel(drep)
        repDir = fullfile(inputDir, drep(i1).name);    
        dtri = dir(fullfile(repDir, 'trial-*-*.mat'));        

        for i2 = 1 : numel(dtri)
            if ~isempty(strfind(dtri(i2).name, '_bad'))
                continue;
            end
            
            trialMat = fullfile(repDir, dtri(i2).name);
            load(trialMat); % gives data
            
            if ~isequal(data.params.name, sent)
                continue;
            end

            % == Determine trial type == %
            if ~isempty(strfind(dtri(i2).name, '-1.mat'))
                trialType = 'N';
            elseif ~isempty(strfind(dtri(i2).name, '-2.mat'))
                trialType = 'R';
            else
                error('Unrecognized trial type in file name: %s', dtri(i2).name);
            end

            % == Load ASR data == %
            trialASRDir = strrep(trialMat, '.mat', '_asr');
            if ~isdir(trialASRDir)
                fprintf(2, 'WARNING: Skipping trial %s due to missing ASR directory\n', trialMat);
                continue;
            end

            julianStdOutFN = fullfile(trialASRDir, 'julian_stdout.txt');
            if ~isfile(julianStdOutFN)
                fprintf(2, 'WARNING: Skipping trial %s due to missing julian_stdout.txt\n', trialMat);
                continue;            
            end

            julianWavFN = fullfile(trialASRDir, 'speech.wav');
            if ~isfile(julianWavFN)
                fprintf(2, 'WARNING: Skipping trial %s due to missing speech.wav\n', trialMat);
                continue;            
            end

            rms_ratio = data.rms(:, 2) ./ data.rms(:, 1);

            t_path = which('parse_asr_out');
            if isempty(t_path)
                addpath(ASR_CODE_PATH);
            end
            asrPAlign = parse_asr_out(julianStdOutFN,  julianWavFN);

            if asrPAlign.nphns ~= nPhones
                fprintf(2, 'WARNING: Skipping trial %s due to erroneous number of phones (possibly speech error)\n', trialMat);
                continue;
            end

            if ~isequal(asrPAlign.phones(2 : end - 1), sentPhns)
                fprintf(2, 'WARNING: Skipping trial %s due to erroneous number of phones (possibly speech error)\n', trialMat);
                continue;
            end

            if ~(isequal(asrPAlign.phones{1}, 'sil') && isequal(asrPAlign.phones{end}, 'sil'))
                fprintf(2, 'WARNING: Skipping trial %s because sil are not found at the beginning and the end of the utterance\n', trialMat);
                continue;
            end

            asrErrCode = check_asrPAlign(asrPAlign);
            if asrErrCode > 0
                fprintf(2, 'WARNING: Skipping trial %s because error detected by check_asrPAlign (errCode = %d)\n', trialMat, asrErrCode);
                continue;
            end

            fprintf(1, 'Processing file: %s (%s)...\n', trialMat, trialType);
            
%             len_s_t_eh_d.(trialType)(end + 1) = asrPAlign.tend(7) - asrPAlign.tbeg(4);
            len_eh_d.(trialType)(end + 1) = asrPAlign.tend(7) - asrPAlign.tbeg(6);

            frameDur = data.params.frameLen / data.params.sr;
            N = size(data.rms, 1);
            tAxis = 0 : frameDur : frameDur * (N - 1);

            idxSpeech = find(tAxis >= asrPAlign.tbeg(2) & tAxis <= asrPAlign.tbeg(end));
            t_rms = data.rms(idxSpeech, 1);

            t_rms_slope = data.rms_slope(idxSpeech);

            if bPlot
                clf;
                for k1 = 1 : 3
                    subplot(3, 1, k1);
                    hold on;

                    if k1 == 1
                        plot(tAxis(idxSpeech), t_rms, 'o-', 'Color', typeColors.(trialType));
                        ylabel('rms');
                    elseif k1 == 2
                        plot(tAxis(idxSpeech), rms_ratio(idxSpeech), 'o-', 'Color', typeColors.(trialType));
                        ylabel('rms\_ratio');
                    else
                        plot(tAxis(idxSpeech), t_rms_slope, 'o-', 'Color', typeColors.(trialType));
                        hold on;
                        plot(tAxis(idxSpeech), zeros(size(idxSpeech)), 'k-');
                        ylabel('rms\_slope');
                    end

                    % == Show phone boundaries == %
                    ys = get(gca, 'YLim');
                    for k2 = 1 : numel(sentPhns)
                        plot(repmat(asrPAlign.tbeg(k2 + 1), 1, 2), ys, 'k-');
                        text(asrPAlign.tbeg(k2 + 1), ys(2) - 0.05 * range(ys), sentPhns{k2}, 'Color', 'k');
                    end
                end
            end

            % === Cumulating statistics for OST parameter calculation === %
            if isequal(sent, 'The steady bat gave birth to pups')
                % == Data for  == %
                steady_s_dur.(trialType)(end + 1) = asrPAlign.tend(4) - asrPAlign.tbeg(4);

                % == For mode 5 (1-2): entering the ah in "the" == %
                idx_dh_ah = find(tAxis >= asrPAlign.tbeg(2) & tAxis < asrPAlign.tend(3));

                rms_dh_ah = data.rms(idx_dh_ah, 1);

                min_rms_dh_ah.(trialType)(end + 1) = min(rms_dh_ah);
                max_rms_dh_ah.(trialType)(end + 1) = max(rms_dh_ah);

                % == For mode 30 (3-4): entering the s in "steady"   %
                %    and mode 31 (5-6): leaving the s in "steady" == %
                idx_s_t = find(tAxis >= asrPAlign.tbeg(4) & tAxis < asrPAlign.tend(5));
                idx_ah = find(tAxis >= asrPAlign.tbeg(3) & tAxis < asrPAlign.tend(3));

                t_rms_ratio = data.rms(:, 2) ./ data.rms(:, 1);
                rms_ratio_s_t = t_rms_ratio(idx_s_t);
                rms_ratio_ah = t_rms_ratio(idx_ah);

                min_rms_ratio_ah.(trialType)(end + 1) = min(rms_ratio_ah);
                max_rms_ratio_s_t.(trialType)(end + 1) = max(rms_ratio_s_t);

                % == For mode 6 (7-8): entering the eh in "steady" == %
                idx_eh = find(tAxis >= asrPAlign.tbeg(6) & tAxis < asrPAlign.tend(6));

                rms_s_t = data.rms(idx_s_t, 1);
                rms_eh = data.rms(idx_eh, 1);

                min_rms_s_t_eh.(trialType)(end + 1) = min([rms_s_t; rms_eh]);
                max_rms_eh.(trialType)(end + 1) = max(rms_eh);

                % == For mode 11 (8-9): leaving the d in "steady" == %            
                idx_eh_d = find(tAxis >= asrPAlign.tbeg(6) & tAxis < asrPAlign.tend(7));
                idx_d_iy = find(tAxis >= asrPAlign.tbeg(7) & tAxis < asrPAlign.tend(8));

                rms_slp_d_iy = data.rms_slope(idx_d_iy);
                [idx_neg_str_beg, idx_neg_str_end] = get_cont_stretches(rms_slp_d_iy < 0);

                if ~isempty(idx_neg_str_beg)
                    neg_str_lens = idx_neg_str_end - idx_neg_str_beg;
                    [max_neg_str_len, idx_max_len] = max(neg_str_lens);
                    max_neg_str_span = sum(rms_slp_d_iy(idx_neg_str_beg(idx_max_len) : idx_neg_str_end(idx_max_len)));
                    neg_str_init = rms_slp_d_iy(idx_neg_str_beg(idx_max_len));

                    max_d_iy_neg_str_len.(trialType)(end + 1) = max_neg_str_len;
                    max_d_iy_neg_str_span.(trialType)(end + 1) = max_neg_str_span;

                    rms_slp_eh_d = data.rms_slope(idx_eh_d);
                    [idx_neg_str_beg, idx_neg_str_end] = get_cont_stretches(rms_slp_eh_d < 0);
                    neg_str_inits = rms_slp_eh_d(idx_neg_str_beg);
                    idx_keep_str = find(neg_str_inits ~= neg_str_init);
                    idx_neg_str_beg = idx_neg_str_beg(idx_keep_str);
                    idx_neg_str_end = idx_neg_str_end(idx_keep_str);

                    if ~isempty(idx_neg_str_beg)
                        neg_str_lens = idx_neg_str_end - idx_neg_str_beg;
                        [prev_max_neg_str_len, prev_idx_max_len] = max(neg_str_lens);
                        prev_max_neg_str_span = sum(rms_slp_eh_d(idx_neg_str_beg(prev_idx_max_len) : idx_neg_str_end(prev_idx_max_len)));

                        max_eh_d_neg_str_len.(trialType)(end + 1) = prev_max_neg_str_len;
                        max_eh_d_neg_str_span.(trialType)(end + 1) = prev_max_neg_str_span;
                    else
                        max_eh_d_neg_str_len.(trialType)(end + 1) = 0;
                        max_eh_d_neg_str_span.(trialType)(end + 1) = 0;
                    end
                end

                % == For mode 10 (10-11): entering the iy in "steady" == %
                [idx_pos_str_beg, idx_pos_str_end] = get_cont_stretches(rms_slp_d_iy > 0);
                if ~isempty(idx_pos_str_beg)
                    pos_str_lens = idx_pos_str_end - idx_pos_str_beg;
                    [max_pos_str_len, idx_max_len] = max(pos_str_lens);

                    max_d_iy_pos_str_len.(trialType)(end + 1) = max_pos_str_len;
                end

                % === Determine the rmsThresh for eh === %
                idx_s2d = find(tAxis >= asrPAlign.tbeg(4) & tAxis < asrPAlign.tend(7));
                rms_s2d = data.rms(idx_s2d, 1);

                min_rms_s2d.(trialType)(end + 1) = min(rms_s2d);
                max_rms_s2d.(trialType)(end + 1) = max(rms_s2d);
            else
                error('The sentence "%s" is currently not supported', sent);
            end
        end
    end

end

%% Determine the maxIOICfg parameers
trialTypes = {'N', 'R'};

for i1 = 1 : numel(trialTypes)
    tt = trialTypes{i1};
    
%     maxInterval_5_10.(tt) = quantile(len_s_t_eh_d.(tt), 0.9) * 1.1;
    maxInterval_5_10.(tt) = quantile(len_eh_d.(tt), 0.9) * 1.1;
end

%% Determine the perturbation parameters %%
if isequal(sent, 'The steady bat gave birth to pups')
    for i1 = 1 : numel(trialTypes)
        tt = trialTypes{i1};
        
        timeWarpCfg.timeWarp_initStat.(tt) = 4;
%         timeWarp_initStat.(tt) = 8;
        timeWarpCfg.timeWarp_tBegin.(tt) = warpOnsetTime; % ARBITRARY!!
        timeWarpCfg.timeWarp_rate1.(tt) = decelWarpRate; % ARBITRARY!!
        timeWarpCfg.timeWarp_dur1.(tt) = mean(steady_s_dur.(tt));        
        timeWarpCfg.timeWarp_durHold.(tt) = timeWarpCfg.timeWarp_dur1.(tt);
        timeWarpCfg.timeWarp_rate2.(tt) = accelWarpRate; % ARBITRARY!!
        
        rmsThresh.(tt) = quantile(min_rms_s2d.(tt), 0.9) + ...
            (quantile(max_rms_s2d.(tt), 0.1) - quantile(min_rms_s2d.(tt), 0.9)) * 0.1;
    end
end

%% Determine the parameter values %%

if isequal(sent, 'The steady bat gave birth to pups')
    nRules = 7;
    
    for i1 = 1 : numel(trialTypes)
        tt = trialTypes{i1};
        
        stat_0_prm1.(tt) = (quantile(min_rms_dh_ah.(tt), 0.9) + quantile(max_rms_dh_ah.(tt), 0.1)) / 2;
        stat_0_prm2.(tt) = 0.02; % ARBITRARY!!
        
        stat_2_prm1.(tt) = (quantile(min_rms_ratio_ah.(tt), 0.9) + quantile(max_rms_ratio_s_t.(tt), 0.1)) / 2;
        stat_2_prm2.(tt) = 0.02; % ARBITRARY!!
        
        stat_4_prm1.(tt) = stat_2_prm1.(tt);
        stat_4_prm2.(tt) = 0.004; % ARBITRARY!!
        
        stat_6_prm1.(tt) = (quantile(min_rms_s_t_eh.(tt), 0.9) + quantile(max_rms_eh.(tt), 0.1)) / 2;
        stat_6_prm2.(tt) = 0.02;
        
        stat_8_prm1.(tt) = round(mean(max_d_iy_neg_str_len.(tt) / 2));
        stat_8_prm2.(tt) = (quantile(max_eh_d_neg_str_span.(tt), 0.1) + quantile(max_d_iy_neg_str_span.(tt), 0.9)) / 2;
                
        stat_10_prm1.(tt) = round(mean(max_d_iy_pos_str_len.(tt) / 2));
        
        ostMat.(tt) = nan(nRules, 4);
        ostMat.(tt)(1, :) = [0, 5, stat_0_prm1.(tt), stat_0_prm2.(tt)];
        ostMat.(tt)(2, :) = [2, 30, stat_2_prm1.(tt), stat_2_prm2.(tt)];
        ostMat.(tt)(3, :) = [4, 31, stat_4_prm1.(tt), stat_4_prm2.(tt)];
        ostMat.(tt)(4, :) = [6, 6, stat_6_prm1.(tt), stat_6_prm2.(tt)];
        ostMat.(tt)(5, :) = [8, 11, stat_8_prm1.(tt), stat_8_prm2.(tt)];
        ostMat.(tt)(6, :) = [10, 10, stat_10_prm1.(tt), NaN];
        ostMat.(tt)(7, :) = [12, 0, NaN, NaN];
        
        % --- Write ostMat to ost file --- %
        if bRunTest
            ost_fn = [mfilename, sprintf('_%s.ost', tt)];
            write_ost_to_file(ostMat.(tt), DEFAULT_RMS_SLOPE_WIN, ost_fn);
            check_file(ost_fn);
        end
    end
else
    error('The sentence "%s" is currently not supported', sent);
end

%% Test shifts
if bRunTest
    for i1 = 1 : numel(drep)
        repDir = fullfile(inputDir, drep(i1).name);    
        dtri = dir(fullfile(repDir, 'trial-*-*.mat'));

        for i2 = 1 : numel(dtri)
            trialMat = fullfile(repDir, dtri(i2).name);
            load(trialMat); % gives data

            if ~isequal(data.params.name, sent)
                continue;
            end

            % == Determine trial type == %
            if ~isempty(strfind(dtri(i2).name, '-1.mat'))
                trialType = 'N';
            elseif ~isempty(strfind(dtri(i2).name, '-2.mat'))
                trialType = 'R';
            else
                error('Unrecognized trial type in file name: %s', dtri(i2).name);
            end

            % == Load ASR data == %
            trialASRDir = strrep(trialMat, '.mat', '_asr');
            if ~isdir(trialASRDir)
                fprintf(2, 'WARNING: Skipping trial %s due to missing ASR directory\n', trialMat);
                continue;
            end

            julianStdOutFN = fullfile(trialASRDir, 'julian_stdout.txt');
            if ~isfile(julianStdOutFN)
                fprintf(2, 'WARNING: Skipping trial %s due to missing julian_stdout.txt\n', trialMat);
                continue;            
            end

            julianWavFN = fullfile(trialASRDir, 'speech.wav');
            if ~isfile(julianWavFN)
                fprintf(2, 'WARNING: Skipping trial %s due to missing speech.wav\n', trialMat);
                continue;            
            end

            rms_ratio = data.rms(:, 2) ./ data.rms(:, 1);

            t_path = which('parse_asr_out');
            if isempty(t_path)
                addpath(ASR_CODE_PATH);
            end
            asrPAlign = parse_asr_out(julianStdOutFN,  julianWavFN);

            if asrPAlign.nphns ~= nPhones
                fprintf(2, 'WARNING: Skipping trial %s due to erroneous number of phones (possibly speech error)\n', trialMat);
                continue;
            end

            if ~(isequal(asrPAlign.phones{1}, 'sil') && isequal(asrPAlign.phones{end}, 'sil'))
                fprintf(2, 'WARNING: Skipping trial %s because sil are not found at the beginning and the end of the utterance\n', trialMat);
                continue;
            end

            fprintf(1, 'Running test on file: %s (%s)...\n', trialMat, trialType);

            ost_fn = [mfilename, sprintf('_%s.ost', trialType)];

            pipCfg_fn = [mfilename, sprintf('_%s.pcf', trialType)];
            format_pcf(basePCF, pipCfg_fn, ...
                       timeWarpCfg.timeWarp_initStat.(tt), ...
                       timeWarpCfg.timeWarp_tBegin.(tt), ...
                       timeWarpCfg.timeWarp_rate1.(tt), ...
                       timeWarpCfg.timeWarp_dur1.(tt), ...
                       timeWarpCfg.timeWarp_durHold.(tt), ...
                       timeWarpCfg.timeWarp_rate2.(tt), ...
                       F1ShiftRatio, FmtShiftStat0, FmtShiftStat1);
            check_file(pipCfg_fn);

            rhy_pert_test_1(trialMat, ost_fn, pipCfg_fn, 1);

            ca;
            drawnow;
        end
    end
end

return