classdef vocalData < ephysData
    properties(SetAccess = private)
        spikeRange = [-8 2]
        callRange = [-4 2]
        baselineRange = [-8 -4]
        smoothingOffset = 0.1
        dT = 1e-3
        minCalls = 7
        minSpikes = 2
        nStd = 2
        consecBins = 15;
        nBoot = 1000
        sortingMetric = 'sortingQuality'
        sortingThreshold = 2
        latencyType = 'pre'
        latencyRange = [-0.5 0.5]
        manualCorrect = 0
        constantBW = 0.1
        preCall = 0.25;
        postCall = 0.25;
        responsiveAlpha = 0.05
        waveformSamples = 32;
        timeWarp = false
        warp_call_length = 0.1
        externalBaseline = false
        
        isolationDistance
        LRatio
        sortingQuality
        manual_correct_info
        onset_or_offset
        daysOld
        batNum
        cellInfo
        expDay
        callSpikes
        callNum
        callLength
        trialFR
        avgFR
        devFR
        avgBaseline
        devBaseline
        latency
        respType
        respValency
        respStrength
        usable
        time
        nCells
        meanFR
        medianISI
        peak2trough
        spikeWidth
        avg_spike
        duration
        misclassified
    end
    properties (SetAccess = public)
        n_age_bins = 7
    end
    properties (Dependent)
        responsiveCells
        youngestAge
        oldestAge
        ageBins
        responsive_cells_by_bat
        weeks_old_bins
        avgLatency
        stdLatency
        bsStdLatencyCI
    end
    methods
        function vd = vocalData(init_or_update,varargin)
            switch init_or_update
                case 'init'
                    if isempty(varargin)
                        call_echo = 'call';
                        onset_or_offset = 'onset';
                    elseif length(varargin) == 1
                        call_echo = varargin{1};
                        onset_or_offset = 'onset';
                    elseif length(varargin) == 2
                        call_echo = varargin{1};
                        onset_or_offset = varargin{2};
                    elseif length(varargin) == 3
                        call_echo = varargin{1};
                        onset_or_offset = varargin{2};
                        externalInput = varargin{3};
                        if strcmp(externalInput.varName,'baseline');
                            vd.externalBaseline = true;
                        end                            
                    end
                    vd.call_echo = call_echo;
                    vd.onset_or_offset = onset_or_offset;
                case 'update'
                    vd = varargin{1};
            end
            dateFormat = 'yyyyMMdd';
            
            time = vd.callRange(1):vd.dT:vd.callRange(2);
            vd.time = inRange(time,[time(1)+vd.smoothingOffset,time(end)-vd.smoothingOffset]);
            
            nBats = length(vd.batNums);
            sortedCells = load([vd.analysisDir 'sortedCells.mat']);
            sortingInfo = load([vd.analysisDir 'sortingInfo.mat']);
            sortingInfo = sortingInfo.sortingInfo ;
            vd.nCells = length(sortedCells.cellList);
            
            if strcmp(vd.call_echo,'play')
                
            else
                
            end
            
            [vd.isolationDistance, vd.LRatio, vd.sortingQuality, vd.avgBaseline,...
                vd.devBaseline, vd.respValency, vd.respStrength, vd.meanFR,...
                vd.medianISI, vd.peak2trough, vd.spikeWidth, vd.duration,...
                vd.daysOld, vd.latency, vd.respType]  = deal(nan(1,vd.nCells));
            [vd.batNum, vd.cellInfo, vd.expDay, vd.callSpikes, vd.callNum,...
                vd.callLength, vd.trialFR, vd.avgFR, vd.devFR] = deal(cell(1,vd.nCells));
                        

            vd.respType = num2cell(vd.respType);
            vd.usable = false(1,vd.nCells);
            vd.avg_spike = nan(vd.nCells,vd.waveformSamples);

            
            cell_k = 1;
            for b = 1:nBats
                cellList = sortedCells.cellList(strcmp(sortedCells.batNumList,vd.batNums{b}));
                for d = 1:length(cellList)
                    vd.batNum{cell_k} = vd.batNums{b};
                    vd.cellInfo{cell_k} = cellList{d};
                    sortingInfo_idx = strcmp({sortingInfo.cellInfo},vd.cellInfo{cell_k}) & strcmp({sortingInfo.batNum},vd.batNum{cell_k});
                    vd.isolationDistance(cell_k) = sortingInfo(sortingInfo_idx).isolationDistance;
                    vd.LRatio(cell_k) = sortingInfo(sortingInfo_idx).LRatio;
                    vd.sortingQuality(cell_k) = sortingInfo(sortingInfo_idx).sortingQuality;
                    vd.expDay{cell_k} = cellList{d}(1:8);
                    [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.callLength{cell_k}, vd.meanFR(cell_k), vd.medianISI(cell_k),  vd.peak2trough(cell_k), vd.spikeWidth(cell_k), vd.avg_spike(cell_k,:), vd.duration(cell_k)] = getSpikes(vd,b,cell_k);
                    vd.daysOld(cell_k) = days(datetime(vd.cellInfo{cell_k}(1:8),'InputFormat',dateFormat) - vd.birthDates{b});
                    if checkUsability(vd,cell_k)
                        vd.usable(cell_k) = true;
                        [vd.avgFR{cell_k}, vd.devFR{cell_k}, vd.trialFR{cell_k}] = gaussKernelVD(vd,cell_k);
                        if ~vd.externalBaseline
                            [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k);       
                        else
                            vd.avgBaseline(cell_k) = externalInput.avgBaseline(cell_k);
                            vd.devBaseline(cell_k) = externalInput.devBaseline(cell_k);
                        end
                        
                        [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                    end
                    cell_k = cell_k + 1;
                end
            end
        end
        function n = numArgumentsFromSubscript(obj,~,~)
            n = numel(obj);
        end
        function varargout = subsref(vd,S)
            if length(S) == 2
                switch S(1).type
                    case '()'
                        nSubs = length(S(1).subs);
                        if ~rem(nSubs,2)
                            cellIdx = true(1,vd.nCells);
                            for idx = 1:2:nSubs
                                switch S(1).subs{idx}
                                    case 'cellInfo'
                                        cellIdx = cellIdx & strcmp(vd.cellInfo,S(1).subs{idx+1});
                                    case 'daysOld'
                                        cellIdx = cellIdx & vd.daysOld==S(1).subs{idx+1};
                                    case 'expDay'
                                        cellIdx = cellIdx & strcmp(vd.expDay,S(1).subs{idx+1});
                                    case 'batNum'
                                        cellIdx = cellIdx & strcmp(vd.batNum,S(1).subs{idx+1});
                                    case 'respType'
                                        cellIdx = cellIdx & strcmp(vd.respType,S(1).subs{idx+1});
                                    case 'respValency'
                                        cellIdx = cellIdx & vd.respValency==S(1).subs{idx+1};
                                    case 'cell_k'
                                        cellIdx = cellIdx & ismember(1:vd.nCells,S(1).subs{idx+1});
                                    otherwise
                                        display('indexing variable not recognized')
                                        return
                                end
                            end
                            if iscell(vd.(S(2).subs)(cellIdx)) && ~any(cellfun(@ischar,vd.(S(2).subs)(cellIdx)))
                                try
                                    varargout = {vertcat(vd.(S(2).subs){cellIdx})};
                                catch
                                    try
                                        varargout = {[vd.(S(2).subs){cellIdx}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                                
                            else
                                varargout = {vd.(S(2).subs)(cellIdx)};
                            end
                        else
                            display('Indexing in VocalData must come in pairs');
                            return
                        end
                    otherwise
                        switch S(2).type
                            case '{}'
                                try
                                    varargout = {vertcat(vd.(S(1).subs){S(2).subs{:}})};
                                catch
                                    try
                                        varargout = {[vd.(S(1).subs){S(2).subs{:}}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                            otherwise
                                try
                                    varargout = {builtin('subsref',vd,S)};
                                catch err
                                    switch err.message
                                        case 'Too many output arguments.'
                                            builtin('subsref',vd,S);
                                        otherwise
                                            display(err)
                                            return
                                    end
                                end
                        end
                end
            else
                try
                    varargout = {builtin('subsref',vd,S)};
                catch err
                    switch err.message
                        case 'Too many output arguments.'
                            builtin('subsref',vd,S);
                        otherwise
                            display(err)
                            return
                    end
                end
            end
            
            
        end
        function vd = updateVar(vd,updateVar,updateVal)
            if ~any(strcmp(updateVar,properties(vd)));
                display('Unrecognized property');
                return
            end
            vd.(updateVar) = updateVal;
            if any(strcmp(updateVar,{'spikeRange','dT','constantBW','callRange','timeWarp','warp_call_length'}))
                vd = vocalData('update',vd);
            elseif any(strcmp(updateVar,{'minCalls','minSpikes','preCall','postCall','sortingMetric','sortingThreshold'}))
                vd = updateUsability(vd);
            elseif strcmp(updateVar,'baselineRange');
                vd = updateBaseline(vd);
            elseif any(strcmp(updateVar,{'nStd','consecBins','latencyType','latencyRange','responsiveAlpha'}))
                vd = updateLatency(vd);
            end
            
        end
        function vd = checkClass(vd,varargin)
            if isempty(varargin)
                vd.misclassified = nan(1,vd.nCells);
                c = vd.responsiveCells(1);
            else
                c = varargin{1};
            end
            figure('units','normalized','outerposition',[0 0 1 1]);
            hold on
            for c = vd.responsiveCells(find(vd.responsiveCells==c):end)
                plotFR(vd,c)
                commandwindow
                class = input('misclassified?');
                if class == 1 || class == 0
                    vd.misclassified(c) = class;
                else
                    vd.misclassified(c) = Inf;
                    break
                end
                cla
            end
            
            
        end
        
        function plotFR(vd,cell_k)
            tRange = [min(vd.time) max(vd.time)];
            if strcmp(vd.call_echo,'call')
                lineColor = 'b';
            elseif strcmp(vd.call_echo,'echo')
                lineColor = 'r';
            end
            boundedline(vd.time,vd.avgFR{cell_k},vd.devFR{cell_k},lineColor);
            if vd.avgBaseline(cell_k) - vd.devBaseline(cell_k) > 0
                plot(tRange,repmat(vd.avgBaseline(cell_k) - vd.devBaseline(cell_k),1,2),'k.-','LineWidth',4)
            else
                plot(tRange,zeros(1,2),'k.-','LineWidth',4)
            end
            plot(tRange,repmat(vd.avgBaseline(cell_k) + vd.devBaseline(cell_k),1,2),'k.-','LineWidth',4)
            plot(repmat(vd.latency(cell_k),1,2),get(gca,'ylim'),'k','LineWidth',5);
            xlim(tRange);
            
            xlabel('Time (s)');
            ylabel('Firing Rate (Hz)');
        end
        function plotRaster(vd,cell_k,varargin)
            hold on
            tRange = [min(vd.time) max(vd.time)];
            if ~isempty(varargin)
                order = varargin{1};
                if size(order,1) ~= 1
                    order = num2cell(order');
                else
                    order = num2cell(order);
                end
            else
                order = num2cell(1:length(vd.callSpikes{cell_k}));
            end
            cellfun(@(trial,r) plot(r,repmat(trial,length(r)),'k.','MarkerSize',8),order,vd.callSpikes{cell_k})
            xlim(tRange);
            ylim([0 length(vd.callSpikes{cell_k})+1]);
            xlabel('time (s)');
            ylabel('trial #');
        end
        function generateAllPlots(vd,saveDir)
            for c = find(vd.usable)
                fName = [saveDir vd.batNum{c} '_' vd.cellInfo{c}];
                h = figure;
                subplot(2,1,1)
                title([vd.batNum{c} '--' vd.cellInfo{c}]);
                vd.plotRaster(c)
                subplot(2,1,2)
                vd.plotFR(c);
                pause(0.00001);
                frame_h = get(h,'JavaFrame');
                set(frame_h,'Maximized',1);
                drawnow;
                saveas(h,fName,'jpg')
                close(h);
            end
        end
        
        function idx = callDataIdx(vd,cData,cell_k)
            firstCall_k = find(strcmp(cData.expDay,vd.cellInfo{cell_k}(1:8)) & strcmp(cData.batNum,vd.batNum{cell_k}),1,'first');
            idx = firstCall_k - 1 + [vd.callNum{cell_k}];
        end
        
        function varargout = scCorr(vd,cell_1,cell_2,period)
            switch period
                case 'pre'
                    spikes1 = cellfun(@(spikes) length(inRange(spikes,[0,-vd.preCall]))/vd.preCall,vd.callSpikes{cell_1});
                    spikes2 = cellfun(@(spikes) length(inRange(spikes,[0,-vd.preCall]))/vd.preCall,vd.callSpikes{cell_2});
                    [spikeCorr, p] = corrcoef(spikes1,spikes2);
                    spikeCorr = spikeCorr(1,2);
                    p = p(1,2);
                case 'post'
                    spikes1 = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_1});
                    spikes2 = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_2});
                    [spikeCorr, p] = corrcoef(spikes1,spikes2);
                    spikeCorr = spikeCorr(1,2);
                    p = p(1,2);
                case 'during'
                    spikes1 = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_1},num2cell(vd.callLength{cell_1}));
                    spikes2 = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_2},num2cell(vd.callLength{cell_2}));
                    [spikeCorr, p] = corrcoef(spikes1,spikes2);
                    spikeCorr = spikeCorr(1,2);
                    p = p(1,2);
                case 'base'
                    spikes1 = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/abs(diff(vd.baselineRange)),vd.callSpikes{cell_1});
                    spikes2 = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/abs(diff(vd.baselineRange)),vd.callSpikes{cell_2});
                    [spikeCorr, p] = corrcoef(spikes1,spikes2);
                    spikeCorr = spikeCorr(1,2);
                    p = p(1,2);
            end
            varargout = {spikeCorr, p};
        end
        
        function callFR = getCallFR(vd,cell_k)
            callFR = cellfun(@(spikes,callLength) length(inRange(spikes, [0,callLength]))/callLength,vd.callSpikes{cell_k},num2cell(vd.callLength{cell_k}));
        end
        
        function responsiveCells = get.responsiveCells(obj)
            responsiveCells = find(~isnan(obj.latency));
        end
        function oldestAge = get.oldestAge(obj)
            oldestAge = max(obj.daysOld);
        end
        function youngestAge = get.youngestAge(obj)
            youngestAge = min(obj.daysOld);
        end
        function ageBins = get.ageBins(obj)
            ageBins = linspace(obj.youngestAge,obj.oldestAge+1,obj.n_age_bins);
        end
        function weeks_old_bins = get.weeks_old_bins(obj)
            weeks_old_bins = obj.ageBins(1:end-1)/7;
        end
        function avgLatency = get.avgLatency(obj)
            avgLatency = zeros(1,length(obj.ageBins)-1);
            for w = 1:length(obj.ageBins)-1
                [~,idx] = inRange(obj.daysOld,[obj.ageBins(w),obj.ageBins(w+1)]);
                latency_in_bin = obj.latency(idx);
                avgLatency(w) = nanmean(latency_in_bin);
            end
        end
        function stdLatency = get.stdLatency(obj)
            stdLatency = zeros(1,length(obj.ageBins)-1);
            for w = 1:length(obj.ageBins)-1
                [~,idx] = inRange(obj.daysOld,[obj.ageBins(w),obj.ageBins(w+1)]);
                latency_in_bin = obj.latency(idx);
                %                stdLatency(w) = nanstd(bootstrp(obj.nBoot,@nanmean,latency_in_bin));
                stdLatency(w) = nanstd(latency_in_bin);
            end
        end
        
        function bsStdLatencyCI = get.bsStdLatencyCI(obj)
            bsStdLatencyCI = zeros(2,length(obj.ageBins)-1);
            for w = 1:length(obj.ageBins)-1
                [~,idx] = inRange(obj.daysOld,[obj.ageBins(w),obj.ageBins(w+1)]);
                latency_in_bin = obj.latency(idx);
                bsStdLatency = bootstrp(obj.nBoot,@nanstd,latency_in_bin);
                bsStdLatencyCI(1,w) = 2*obj.stdLatency(w) - quantile(bsStdLatency,0.975);
                bsStdLatencyCI(2,w) = 2*obj.stdLatency(w) - quantile(bsStdLatency,0.025);
            end
            
        end
        function vd = set_manual_responsive(vd)
            
            try
                s = load([vd.analysisDir 'manual_responsive_cell_correction_' vd.call_echo '.mat']);
                [e1{1:length(s.manualCorrectFN)}] = deal('falseNegative');
                [e2{1:length(s.manualCorrectFP)}] = deal('falsePositive');
                manualCorrect_cellInfo = [s.manualCorrectFN s.manualCorrectFP];
                manualCorrect_batNum = [s.manualCorrectFN_batnum s.manualCorrectFP_batnum];
                errorTypes = [e2 e1];
                cell_k = cell(1,length(errorTypes));
                for c = 1:length(errorTypes)
                    errorType = errorTypes{c};
                    mc_cellInfo = manualCorrect_cellInfo{c};
                    mc_batNum = manualCorrect_batNum{c};
                    cell_k{c} = find(strcmp(vd.cellInfo,mc_cellInfo) & strcmp(vd.batNum,mc_batNum));
                    vd = setResponsive(vd,errorType,cell_k{c});
                end
            catch
                display('couldn not find manual correction data')
                for c = vd.responsiveCells
                    
                end
            end
            
            
            manual_correct = struct('cell_k',cell_k,'errorType',errorTypes,'latency',num2cell([nan(1,length(s.manualCorrectFP)) vd.latency([cell_k{:}])]));
            save([vd.analysisDir 'manual_correct_' vd.call_echo '.mat'],'manual_correct');
            vd = load_manual_responsive(vd,'load');
            
            
        end
        function vd = load_manual_responsive(vd,loadFlag)
            
            switch loadFlag
                case 'unload'
                    vd.manualCorrect = 0;
                    vd.manual_correct_info = [];
                    vd = updateLatency(vd);
                case 'load'
                    if strcmp(vd.call_echo,'call')
                        s = load([vd.analysisDir 'manual_correct_call.mat']);
                    else
                        s = load([vd.analysisDir 'manual_correct_echo.mat']);
                    end
                    vd.manualCorrect = 1;
                    vd.manual_correct_info = s.manual_correct;
                    vd = updateManual(vd);
            end
            
        end
        
        function responsive_cells_by_bat = get.responsive_cells_by_bat(obj)
            responsive_cells_by_bat = cell(1,length(obj.batNums));
            for n = 1:length(obj.batNums)
                responsive_cells_by_bat{n} = intersect(find(strcmp(obj.batNums{n},obj.batNum)),obj.responsiveCells);
            end
        end
        function lm = fitLatencyAgeLM(obj)
            lm = fitlm(obj.daysOld(obj.responsiveCells),obj.latency(obj.responsiveCells));
        end
        
    end
    
end

function [callSpikes, callNum, callLength, meanFR, medianISI, peak2trough, spikeWidth, avg_spike, duration] = getSpikes(vd,b,cell_k)

cellInfo = vd.cellInfo{cell_k};
baseDir = vd.baseDirs{b};
batNum = vd.batNums{b};
expDate = cellInfo(1:8);
tetrodeNum = str2double(cellInfo(11));
cellNum = str2double(cellInfo(16:17));
call_echo = vd.call_echo;
spikeRange = 1e3*vd.spikeRange;
call_onset_offset = vd.onset_or_offset;

if cellNum < 10
    sorted_cell_string = '_SS_0';
else
    sorted_cell_string = '_SS_';
end

tetrode = [expDate 'TT' num2str(tetrodeNum) sorted_cell_string num2str(cellNum)]; % build tetrode/cell string
ttDir = ['C:\Users\phyllo\Documents\Maimon\ephys\Data_processed\bat' batNum filesep expDate filesep tetrode '.ntt']; % filename for cell of interest
audioDir = [baseDir 'bat' batNum filesep 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
load('C:\Users\phyllo\Documents\Maimon\ephys\data_analysis_results\cell_stability_info');
idx = find(strcmp(tetrode,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
timeWarp = vd.timeWarp;
warp_call_length = vd.warp_call_length;


if strcmp(call_echo,'call') % are we looking at calls or echolocation clicks?
    callStr = 'calls';
    analyzed_audio_dir = [audioDir 'Analyzed_auto\']; % point to pre-cut calls
    s = load([analyzed_audio_dir 'juv_calls.mat']);
    juv_calls = s.juv_calls; % extract data from loaded structure
    callFiles = juv_calls ~= -1; % where juv_call == -1 is to note that the cutting procedure erroneously included a sound that was not a call
elseif strcmp(call_echo,'echo')
    callStr = 'echoes';
    analyzed_audio_dir = [audioDir 'Analyzed_auto_echo\']; % point to pre-cut echo. clicks
    callFiles = true(1,1000); % so initialize echolocation files 1000 total clicks (will be truncated later, but could cause issues if total number of clicks > 1000)
end

audio2nlg = load([audioDir 'audio2nlg_fit.mat']); % load fit data to sync audio to nlg data
stabilityBounds = stabilityBounds - audio2nlg.first_nlg_pulse_time;
timestamps = Nlx2MatSpike(ttDir,[1 0 0 0 0],0,1,[]); % load sorted cell data
timestamps = 1e-3*timestamps - audio2nlg.first_nlg_pulse_time; % convert to ms and align to first TTL pulse on the NLG
timestamps = inRange(timestamps,stabilityBounds);

[callSpikes, callNum, callLength] = get_call_spikes(timestamps,stabilityBounds,audioDir,analyzed_audio_dir,callFiles,callStr,spikeRange,call_onset_offset,timeWarp,warp_call_length);
[meanFR, medianISI, peak2trough, spikeWidth, avg_spike, duration] = get_spike_stats(timestamps,ttDir);
end

function [callSpikes, callNum, callLength] = get_call_spikes(timestamps,stabilityBounds,audioDir,analyzed_audio_dir,callFiles,callStr,spikeRange,call_onset_offset,timeWarp,warp_call_length)

nCalls = length(callFiles); % total # of calls, excluding any errors by cutting procedure

% get_corrected_call_times:
% function to return: 1) a structure containing a) the number of the .wav file that each call was taken from and b) call times synched to the NLG
% and 2) all synced calls times
[cut_call_data, all_call_times] = get_corrected_call_times(audioDir,analyzed_audio_dir,callStr);
if (length(callFiles) ~= length(all_call_times)) && strcmp(callStr,'echoes') && ~isempty(all_call_times)
    callFiles = callFiles(1:length(all_call_times)/2);
elseif isempty(all_call_times)
    callFiles = [];
end
all_call_times = all_call_times(repmat(callFiles, [1 2])); % remove non-call files from list of call times

n_cut_call_files = length(cut_call_data); % total number of .wav files with calls

call_k = 1; % counter for all used calls
rec_k = 1; % counter for all calls (to exclude any non-call recordings)
rec_file_k = 1; % counter for all call files (including non-call recordings)

callSpikes = cell(1,nCalls); % initialize cell to contain spike times relative to call onset
callNum = nan(1,nCalls);
callLength = nan(1,nCalls);

for call_f = 1:n_cut_call_files % itereate through .wav files with calls
    n_call_f = length(cut_call_data(call_f).cutcalls);
    for call = 1:n_call_f % iterate through all the calls within this .wav file
        cp = cut_call_data(call_f).corrected_callpos(call,:); % time of onset and offset of the call of interest on this iteration
        callposOffset = [cp(1)+spikeRange(1), cp(2)+spikeRange(2)]; % total time window around call of interest
        
        if any(rec_file_k == find(~callFiles)) || any(isnan(cut_call_data(call_f).corrected_callpos(call,:))) % skip non-call recordings or recordings before TTL alignment
            rec_file_k = rec_file_k + 1;
            continue
        elseif ~(callposOffset(1) > stabilityBounds(1) && callposOffset(2) < stabilityBounds(2)) ||... % check that this call is within timeframe that this cell is stable
                (strcmp(call_onset_offset,'onset') && any(all_call_times > (cp(1) + spikeRange(1)) & all_call_times < cp(1))) || ... skip any calls within a call train
                (strcmp(call_onset_offset,'offset') && any(all_call_times < (cp(2) + spikeRange(2)) & all_call_times > cp(2)))
            rec_k = rec_k + 1;
            rec_file_k = rec_file_k + 1;
            continue
        end
        
        relativeCP = cp(strcmp({'onset','offset'},call_onset_offset));
        
        timestampsCall = inRange(timestamps,callposOffset); % all spikes occuring within that time window
        callSpikes{call_k} = 1e-3*(timestampsCall - relativeCP); % align spikes to call onset, convert to sec, and store
        callNum(call_k) = rec_k;
        callLength(call_k) = 1e-3*diff(cp);
        if timeWarp
            callSpikes{call_k} = callSpikes{call_k} * warp_call_length/callLength(call_k);
        end
        call_k = call_k + 1; % increment total calls
        rec_k = rec_k + 1; % increment total recordings
        rec_file_k = rec_file_k + 1;
    end
end

%now truncate the list of calls/clicks to how many calls/clicks actually occurred
nCalls = call_k - 1; % determine the actual number of calls
callSpikes = callSpikes(1:nCalls); % truncate list of spikes relative to call onset to length of n_calls
callNum = callNum(1:nCalls);
callLength = callLength(1:nCalls);
end

function [baselineMean, baselineStd] = calculate_baseline(vd,cell_k)

baselineLength = abs(diff(vd.baselineRange));
nTrial = length(vd.callSpikes{cell_k});
allSpikes = [vd.callSpikes{cell_k}{:}];
usedSpikes = inRange(allSpikes,vd.baselineRange);
baselineMean = length(usedSpikes)/nTrial/baselineLength;
baselineStd = std(cellfun(@(x) length(inRange(x,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}));

end

function usable = checkUsability(vd,cell_k)

% usable = all(cellfun(@(x) length(x(x > vd.callRange(1) & x < vd.callRange(2))),vd.callSpikes{cell_k})>=vd.minSpikes) && length(vd.callSpikes{cell_k}) >= vd.minCalls;
allSpikes = vd.callSpikes{cell_k};
allSpikes = [allSpikes{:}];
preSpikes = inRange(allSpikes,[-vd.preCall,0]);
postSpikes = inRange(allSpikes,[0,vd.postCall]);
nTrial = size(vd.callSpikes{cell_k},2);
switch vd.sortingMetric
    case 'sortingQuality'
        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
    case 'ratio_and_distance'
        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
end
usable =  length(preSpikes) >= vd.minSpikes &&...
    length(postSpikes) >= vd.minSpikes &&...
    length(vd.callSpikes{cell_k}) >= vd.minCalls &&...
    wellSorted;

end

function [spike_rate_mean, spike_rate_sem, spikeRate, bw] = gaussKernelVD(vd,cell_k)

time = vd.callRange(1):vd.dT:vd.callRange(2);

nTrial = length(vd.callSpikes{cell_k});
nT = length(time);
if isnan(vd.constantBW)
    bw = kde([vd.callSpikes{cell_k}{:}]);
else
    bw = vd.constantBW;
end
spikeRate = zeros(nTrial,nT);
for tt = 1:nTrial
    kSpike = 1;
    usedSpikes = [vd.callSpikes{cell_k}{tt}];
    usedSpikes = inRange(usedSpikes,[time(1),time(end)]);
    if ~isempty(usedSpikes)
        n_used_spikes = length(usedSpikes);
        spike_rate_trial = zeros(n_used_spikes,nT);
        for s = 1:n_used_spikes
            kernel = Gauss(time,usedSpikes(s),bw.^2);
            spike_rate_trial(kSpike,:) = kernel*vd.dT;
            kSpike = kSpike + 1;
        end
        spikeRate(tt,:) = sum(spike_rate_trial,1)/vd.dT;
    end
    
end

[~,idx] = inRange(time,[time(1)+vd.smoothingOffset,time(end)-vd.smoothingOffset]);
spikeRate = spikeRate(:,idx);
spike_rate_mean = mean(spikeRate);
spike_rate_sem = smooth(std(spikeRate)/sqrt(nTrial))';

end

function y = Gauss(t,tS,w)

y = (1./sqrt(2*pi*w)) * exp( -(((t - tS).^2)./(2*w)));

end

function [latency, respType, respValency, respStrength] = calculateLatency(vd,cell_k)

[responsive, respType, respValency, respStrength] = periCallResponsive(vd,cell_k);

if responsive
    switch respType
        case 'pre'
            [~,idx1] = inRange(vd.time,[0 -Inf]);
            [~,idx2] = inRange(vd.time,vd.latencyRange);
            idx = idx1 & idx2;
            latencyTime = vd.time(idx);
            latencyFR = vd.avgFR{cell_k}(idx);
        case 'post'
            [~,idx1] = inRange(vd.time,[0,Inf]);
            [~,idx2] = inRange(vd.time,vd.latencyRange);
            idx = idx1 & idx2;
            latencyTime = vd.time(idx);
            latencyFR = vd.avgFR{cell_k}(idx);
        case 'during'
            [~,idx] = inRange(vd.time, [0, mean([vd.callLength{cell_k}])]);
            latencyTime = vd.time(idx);
            latencyFR = vd.avgFR{cell_k}(idx);
    end
            
    
    [~, max_fr_idx] = max(abs(latencyFR - vd.avgBaseline(cell_k)));
    latency = latencyTime(max_fr_idx);
else
    latency = nan;
end
end

function [responsive, respType, respValency, respStrength] = periCallResponsive(vd,cell_k)

baselineLength = diff(vd.baselineRange);
baseSpikes = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k});

switch vd.latencyType
    case 'during'
        callSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_k},num2cell(vd.callLength{cell_k}));
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'during';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'pre'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k});
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'pre';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'post'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_k});
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'post';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'prePost'
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k});
        postCallSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_k});
        responsivePre = ttest(baseSpikes,preCallSpikes,'Alpha',vd.responsiveAlpha);
        responsivePost = ttest(baseSpikes,postCallSpikes,'Alpha',vd.responsiveAlpha);
        responsive = responsivePre | responsivePost;
        if responsive
            allSpikes = {preCallSpikes,postCallSpikes};
            responsivity = [responsivePre,responsivePost];
            responseTypes = {'pre','post'};
            respStrength = [abs(mean(preCallSpikes) - mean(baseSpikes)),...
                abs(mean(postCallSpikes) - mean(baseSpikes))];
            
            respStrength(~responsivity) = -Inf;
            [respStrength, prePost] =  max(respStrength);
            respValency = sign(mean(allSpikes{prePost}) - mean(baseSpikes));
            respType = responseTypes{prePost};
        end
        
    case 'peri'
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k});
        postCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[callLength,callLength + vd.postCall]))/vd.postCall,vd.callSpikes{cell_k},num2cell(vd.callLength{cell_k}));
        duringCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_k},num2cell(vd.callLength{cell_k}));
        responsivePre = ttest(baseSpikes,preCallSpikes,'Alpha',vd.responsiveAlpha);
        responsivePost = ttest(baseSpikes,postCallSpikes,'Alpha',vd.responsiveAlpha);
        responsiveDuring = ttest(baseSpikes,duringCallSpikes,'Alpha',vd.responsiveAlpha);
        responsive = responsivePre | responsivePost | responsiveDuring;
        if responsive
            allSpikes = {preCallSpikes,postCallSpikes,duringCallSpikes};
            responsivity = [responsivePre,responsivePost,responsiveDuring];
            responseTypes = {'pre','post','during'};
            respStrength = [abs(mean(preCallSpikes) - mean(baseSpikes)),...
                abs(mean(postCallSpikes) - mean(baseSpikes)),...
                abs(mean(duringCallSpikes) - mean(baseSpikes))];
            
            respStrength(~responsivity) = -Inf;
            [respStrength, prePostDuring] =  max(respStrength);
            respValency = sign(mean(allSpikes{prePostDuring}) - mean(baseSpikes));
            respType = responseTypes{prePostDuring};
        end
end

if ~responsive
    respType = NaN;
    respValency = NaN;
    respStrength = NaN;
end


end

function [meanFR, medianISI, peak2trough, spikeWidth, avg_spike, duration] = get_spike_stats(timestamps,ttDir)
ADC_period = 34.13; % in microseconds

[~, samples] = Nlx2MatSpike(ttDir,[1 0 0 0 1],0,1,[]);

dt_height = 0.01;
n_samples = size(samples,1); 
t_height = 1:dt_height:n_samples;

meanFR = length(timestamps)/(1e-3*range(timestamps));
ISIs = diff(timestamps);

% find channel with largest amplitude
[~,spikeChannel] = max(max(mean(samples,3),[],1));
avg_spike_max_tt = mean(squeeze(samples(:,spikeChannel,:)),2);
% set time through waveform and interpolate points to make upsampled waveform 
time = ADC_period*(0:length(avg_spike_max_tt)-1);
interp_time = 0:1:(ADC_period*(length(avg_spike_max_tt)-1));
% upsample waveform
avg_spike_interp = interp1(time,avg_spike_max_tt,interp_time,'spline');
% find peak and trough values and indices
[spikePeak, peakIdx] = max(avg_spike_interp);
[spikeTrough, troughIdx] = min(avg_spike_interp(1:peakIdx));
% find avg waveform on most active channel
avg_spike = avg_spike_max_tt;
% spikeHeight = (max(avg_spike_max_tt) - min(avg_spike_max_tt(1:find(avg_spike_max_tt==max(avg_spike_max_tt)))))/2;
halfSpikePeak = spikePeak/2;
spikeWidth = diff(polyxpoly(1:n_samples,avg_spike_max_tt,t_height,repmat(halfSpikePeak,size(t_height))));
if isempty(spikeWidth)
    spikeWidth = nan;
elseif length(spikeWidth) > 1
    spikeWidth = spikeWidth(1);
end
if troughIdx > 1
    peak2trough = abs(spikeTrough/spikePeak);
    duration = interp_time(troughIdx) - interp_time(peakIdx);
else
    peak2trough = nan;
    duration = nan;
end

medianISI = median(ISIs);
end

function [xSub,idx] = inRange(x,bounds)

bounds = sort(bounds);
idx = x>bounds(1) & x<= bounds(2);
xSub = x(idx);


end

function vd = setResponsive(vd,errorType,cell_k)


if isempty(cell_k)
    display('couldn''t find specified cell')
    keyboard
    return
end
if strcmp(errorType,'falseNegative') && ~ismember(cell_k,vd.responsiveCells)
    figure;
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    title([vd.batNum{cell_k} '--' vd.cellInfo{cell_k}]);
    vd.plotRaster(cell_k)
    subplot(2,1,2)
    vd.plotFR(cell_k);
    %     pause(0.00001);
    %     frame_h = get(h,'JavaFrame');
    %     set(frame_h,'Maximized',1);
    %     drawnow;
    [latencyManual,~] = ginput(1);
    vd.latency(cell_k) = latencyManual;
    cla;
    vd.plotFR(cell_k);
    pause;
    close(h);
elseif strcmp(errorType,'falsePositive') && ismember(cell_k,vd.responsiveCells)
    vd.latency(cell_k) = nan;
else
    display('error type and response profile mismatch')
    keyboard;
end
end

function vd = updateManual(vd)
for c = 1:length(vd.manual_correct_info)
    if ~isfield(vd.manual_correct_info,'cellInfo') && isfield(vd.manual_correct_info,'cell_k')
        cell_k = vd.manual_correct_info(c).cell_k;
    elseif isfield(vd.manual_correct_info,'cellInfo') && isfield(vd.manual_correct_info,'batNum')
        cell_k = find(strcmp(vd.cellInfo,vd.manual_correct_info(c).cellInfo) & strcmp(vd.batNum,vd.manual_correct_info(c).batNum));
    else
        display('couldn''t find appropriate manual correct info file');
        return
    end
    vd.latency(cell_k) = vd.manual_correct_info(c).latency;
end
vd.latency(vd.misclassified==1) = NaN;
end

function vd = updateLatency(vd)
vd.latency = nan(1,vd.nCells);
vd.respType = num2cell(nan(1,vd.nCells));
vd.respValency = nan(1,vd.nCells);
for cell_k = 1:vd.nCells
    if vd.usable(cell_k)
        [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
    end
end
if vd.manualCorrect
    vd = updateManual(vd);
end
end

function vd = updateUsability(vd)

updated_cells = 0;
for cell_k = 1:vd.nCells
    usability = checkUsability(vd,cell_k);
    if vd.usable(cell_k) ~= usability
        if usability - vd.usable(cell_k) == 1 % cell is now usable
            [vd.avgFR{cell_k}, vd.devFR{cell_k}] = gaussKernelVD(vd,cell_k);
            [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k);
        elseif usability - vd.usable(cell_k) == -1 % cell is now unusable
            vd.avgFR{cell_k} = [];
            vd.devFR{cell_k} = [];
            vd.avgBaseline(cell_k) = nan;
            vd.devBaseline(cell_k) = nan;
        end
        vd.usable(cell_k) = usability;
        updated_cells = updated_cells + 1;
    end
end
vd = updateLatency(vd);
display([num2str(updated_cells) ' cells updated']);
end

function vd = updateBaseline(vd)

for cell_k = 1:vd.nCells
    if vd.usable(cell_k)
        [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k);
    end
end
vd = updateLatency(vd);
end
