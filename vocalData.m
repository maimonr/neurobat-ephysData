classdef vocalData < ephysData
    properties(SetAccess = public)
        spikeRange = [-2 2]
        callRange = [-2 2]
        baselineRange = [-2 -1]
        smoothingOffset = 0
        frBandwidthRange = 2e-3:5e-3:1e-1
        kernelType = 'manual'
        dT = 5e-3
        minCalls = 10
        minSpikes = 2
        nStd = 2
        consecBins = 15;
        nBoot = 1000
        sortingMetric = 'sortingQuality'
        sortingThreshold = 2
        latencyType = 'std_above_baseline'
        latencyRange = [-1 1]
        manualCorrect = 0
        constantBW = 0.05
        preCall = 0.5
        postCall = 0.5
        responsiveAlpha = 0.05
        responsive_nStd_over_baseline = 2
        min_responsive_trials = 10
        min_fraction_responsive_trials = 0.25
        min_responsive_length = 0.05
        cusumDelta = 2
        cusum_nStd = 20
        cusumBaseline = -0.05
        cusumStart = -0.05
        waveformSamples = 32;
        timeWarp = false
        warp_call_length = 0.1
        baselineMethod = 'randomSamples'
        n_baseline_samples = 1e2
        baseline_sample_length = 5
        tetrodeStr = 'TT';
        clusterStr = '_SS_';
        
        usedCalls
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
        frBandwidth
        trialFR
        trial_spike_train
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
        tetrodeNum
        tetrodeDepth
    end
    properties (SetAccess = public)
        n_age_bins = 7
        manual_responsive_cells = []
        manual_latency
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
        function vd = vocalData(varargin)
            
            call_echo = 'call';
            onset_or_offset = 'onset';
            expType = 'juvenile';
            updateFlag = false;
            use_select_cells = false;
            
            for n = 1:2:nargin
                switch varargin{n}
                    
                    case 'update'
                        vd_update = varargin{n+1};
                        call_echo = vd_update.call_echo;
                        onset_or_offset = vd_update.onset_or_offset;
                        expType = vd_update.expType;
                        updateFlag = true;
                    
                    case 'call_echo'
                        call_echo = varargin{n+1};
                        
                    case 'onset_or_offset'
                        onset_or_offset = varargin{n+1};
                                            
                    case 'selectCells'
                        use_select_cells = true;
                        selectCells = varargin{n+1};
                        
                    case 'expType'
                        expType = varargin{n+1};
                        
                end
                
            end
            
            vd = vd@ephysData(expType);
            
            if updateFlag
               vd = vd_update; 
            end
            
            vd.call_echo = call_echo;
            vd.onset_or_offset = onset_or_offset;
            
            time = vd.callRange(1):vd.dT:vd.callRange(2);
            vd.time = inRange(time,[time(1)+vd.smoothingOffset,time(end)-vd.smoothingOffset]);
            
            vd.expDay = datetime([],[],[]);            
            
            cell_k = 1;
            lastProgress = 0;

            switch vd.expType
                
                case 'adult'
                    
                    spike_file_names = dir([vd.spike_data_dir '*.csv']);
                    spike_file_names = spike_file_names(arrayfun(@(x) contains(x.name,vd.batNums),spike_file_names));
                    spike_file_name_dlm = '_';
                    spike_file_name_format = 'bbbbb_yyyymmddTTt_SS_ss';
                    nBats = length(vd.batNums);
                    vd.nCells = length(spike_file_names);
                    
                    [vd.avgBaseline, vd.devBaseline, vd.respValency,...
                        vd.respStrength, vd.meanFR, vd.latency, vd.respType,...
                        vd.tetrodeNum, vd.tetrodeDepth]  = deal(nan(1,vd.nCells));
                    [vd.batNum, vd.cellInfo, vd.callSpikes, vd.callNum, vd.frBandwidth,...
                        vd.callLength, vd.trialFR, vd.avgFR, vd.devFR,...
                        vd.trial_spike_train, vd.usedCalls] = deal(cell(1,vd.nCells));
                    vd.respType = num2cell(vd.respType);
                    vd.usable = false(1,vd.nCells);
                    for b = 1:nBats
                        spike_file_names = dir([vd.spike_data_dir '*.csv']);
                        spike_file_names = spike_file_names(arrayfun(@(x) contains(x.name,vd.batNums{b}),spike_file_names));
                        for d = 1:length(spike_file_names)
                            
                            vd.batNum{cell_k} = vd.batNums{b};
                            idx = strfind(spike_file_names(d).name,[vd.batNums{b} spike_file_name_dlm]);
                            cellInfo = spike_file_names(d).name(idx+length([vd.batNums{b} spike_file_name_dlm]):idx+length(spike_file_name_format)-1);
                            if use_select_cells && ~any(strcmp(cellInfo,selectCells))
                                cell_k = cell_k + 1;
                                continue
                            end
                            
                            vd.cellInfo{cell_k} = cellInfo;
                            vd.expDay(cell_k) = datetime(vd.cellInfo{cell_k}(1:8),'InputFormat',vd.dateFormat);
                            [stabilityBounds, cut_call_data, ~, ~, success] = getCellInfo(vd,b,cell_k);
                            timestamps = csvread([vd.spike_data_dir spike_file_names(d).name]);
                            
                            if size(timestamps,1) ~= 1
                               timestamps = timestamps'; 
                            end

                            [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.callLength{cell_k}] ...
                                = get_call_spikes(timestamps,stabilityBounds,cut_call_data,...
                                1e3*vd.spikeRange,vd.onset_or_offset,vd.timeWarp,vd.warp_call_length);
                            
                            vd.usedCalls{cell_k} = get_used_calls(stabilityBounds,cut_call_data,1e3*vd.spikeRange,vd.onset_or_offset);
                            vd.sortingQuality(cell_k) = vd.sortingThreshold;
                            if checkUsability(vd,cell_k) && success
                                vd.usable(cell_k) = true;
                                [vd.avgFR{cell_k}, vd.devFR{cell_k}, vd.trialFR{cell_k},...
                                    vd.trial_spike_train{cell_k}, vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
                                [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k,timestamps,cut_call_data);
                                [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                            end
                            
                            progress = 100*(cell_k/vd.nCells);
                            
                            if mod(progress,10) < mod(lastProgress,10)
                                fprintf('%d %% of cells processed\n',round(progress));
                            end
                            
                            lastProgress = progress;
                            
                            cell_k = cell_k + 1;
                        end
                    end
                    
                
                case 'juvenile'
                    
                    nBats = length(vd.batNums);
                    sortedCells = load([vd.analysisDir 'sortedCells.mat']);
                    sortingInfo = load([vd.analysisDir 'sortingInfo.mat']);
                    sortingInfo = sortingInfo.sortingInfo ;
                    vd.nCells = length(sortedCells.cellList);
                    
                    [vd.isolationDistance, vd.LRatio, vd.sortingQuality, vd.avgBaseline,...
                        vd.devBaseline, vd.respValency, vd.respStrength, vd.meanFR,...
                        vd.medianISI, vd.peak2trough, vd.spikeWidth, vd.duration,...
                        vd.daysOld, vd.latency, vd.respType, vd.tetrodeNum, vd.tetrodeDepth]  = deal(nan(1,vd.nCells));
                    [vd.batNum, vd.cellInfo, vd.callSpikes, vd.callNum, vd.frBandwidth,...
                        vd.callLength, vd.trialFR, vd.avgFR, vd.devFR,...
                        vd.trial_spike_train, vd.usedCalls] = deal(cell(1,vd.nCells));
                    vd.respType = num2cell(vd.respType);
                    vd.usable = false(1,vd.nCells);
                    vd.avg_spike = nan(vd.nCells,vd.waveformSamples);
                    
                    tt_depth_format = '%{MM/dd/yy}D %f %f %f %f';
                    ttString = 'TT';
                    
                    for b = 1:nBats
                        cellList = sortedCells.cellList(strcmp(sortedCells.batNumList,vd.batNums{b}));
                        tetrodeDepths = readtable([vd.baseDirs{b} 'bat' vd.batNums{b} filesep 'tetrode_depths_' vd.batNums{b} '.csv'],'format',tt_depth_format);
                        
                        for d = 1:length(cellList)
                            if use_select_cells && ~any(strcmp(cellList{d},selectCells))
                                cell_k = cell_k + 1;
                                continue
                            end
                            vd.batNum{cell_k} = vd.batNums{b};
                            vd.cellInfo{cell_k} = cellList{d};
                            vd.tetrodeNum(cell_k) = str2double(vd.cellInfo{cell_k}(strfind(vd.cellInfo{cell_k},ttString)+length(ttString)));
                            sortingInfo_idx = strcmp({sortingInfo.cellInfo},vd.cellInfo{cell_k}) & strcmp({sortingInfo.batNum},vd.batNum{cell_k});
                            vd.isolationDistance(cell_k) = sortingInfo(sortingInfo_idx).isolationDistance;
                            vd.LRatio(cell_k) = sortingInfo(sortingInfo_idx).LRatio;
                            vd.sortingQuality(cell_k) = sortingInfo(sortingInfo_idx).sortingQuality;
                            vd.expDay(cell_k) = datetime(vd.cellInfo{cell_k}(1:8),'InputFormat',vd.dateFormat);
                            vd.tetrodeDepth(cell_k) = tetrodeDepths{tetrodeDepths.Date==vd.expDay(cell_k),vd.tetrodeNum(cell_k)+1};
                            
                            [stabilityBounds, cut_call_data, ~, ttDir] = getCellInfo(vd,b,cell_k);
                            timestamps = getSpikes(vd,b,cell_k);
                            
                            [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.callLength{cell_k}] ...
                                = get_call_spikes(timestamps,stabilityBounds,cut_call_data,...
                                1e3*vd.spikeRange,vd.onset_or_offset,vd.timeWarp,vd.warp_call_length);
                            
                            vd.usedCalls{cell_k} = get_used_calls(stabilityBounds,cut_call_data,1e3*vd.spikeRange,vd.onset_or_offset);
                            
                            [vd.meanFR(cell_k), vd.medianISI(cell_k), vd.peak2trough(cell_k),...
                                vd.spikeWidth(cell_k), vd.avg_spike(cell_k,:), vd.duration(cell_k)] = get_spike_stats(timestamps,ttDir);
                            
                            vd.daysOld(cell_k) = days(vd.expDay(cell_k) - vd.birthDates{b});
                            
                            if checkUsability(vd,cell_k)
                                vd.usable(cell_k) = true;
                                [vd.avgFR{cell_k}, vd.devFR{cell_k}, vd.trialFR{cell_k},...
                                    vd.trial_spike_train{cell_k}, vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
                                [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k,timestamps,cut_call_data);
                                [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                            end
                            
                            progress = 100*(cell_k/vd.nCells);
                            
                            if mod(progress,10) < mod(lastProgress,10)
                                fprintf('%d %% of cells processed\n',round(progress));
                            end
                            
                            lastProgress = progress;
                            
                            cell_k = cell_k + 1;
                        end
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
                                        daysOldIdx = false(1,vd.nCells);
                                        for d = S(1).subs{idx+1}
                                            daysOldIdx = daysOldIdx | vd.daysOld==d;
                                        end
                                        cellIdx = cellIdx & daysOldIdx;
                                    case 'expDay'
                                        if length(S(1).subs{idx+1}) == 1
                                            cellIdx = cellIdx & (vd.expDay == S(1).subs{idx+1});
                                        else
                                            cellIdx = cellIdx & (vd.expDay >= S(1).subs{idx+1}(1) & vd.expDay < S(1).subs{idx+1}(2));
                                        end
                                    case 'batNum'
                                        cellIdx = cellIdx & strcmp(vd.batNum,S(1).subs{idx+1});
                                    case 'respType'
                                        cellIdx = cellIdx & strcmp(vd.respType,S(1).subs{idx+1});
                                    case 'respValency'
                                        cellIdx = cellIdx & vd.respValency==S(1).subs{idx+1};
                                    case 'cell_k'
                                        cellIdx = cellIdx & ismember(1:vd.nCells,S(1).subs{idx+1});
                                    case 'misclassified'
                                        cellIdx = cellIdx & vd.misclassified==S(1).subs{idx+1};
                                    case 'responsive'
                                        cellIdx = cellIdx & ~isnan(vd.latency)==S(1).subs{idx+1};
                                    case 'tetrodeNum'
                                        cellIdx = cellIdx & vd.tetrodeNum==S(1).subs{idx+1};
                                    case 'usable'
                                        cellIdx = cellIdx & vd.usable==S(1).subs{idx+1};
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
            if ~any(strcmp(updateVar,properties(vd)))
                disp('Unrecognized property');
                return
            end
            vd.(updateVar) = updateVal;
            if any(strcmp(updateVar,{'spikeRange','dT','constantBW','callRange','timeWarp','warp_call_length','onset_or_offset','kernelType','frBandwidthRange'}))
                vd = vocalData('update',vd);
            elseif any(strcmp(updateVar,{'minCalls','minSpikes','preCall','postCall','sortingMetric','sortingThreshold'}))
                vd = updateUsability(vd);
            elseif any(strcmp(updateVar,{'baselineMethod','n_baseline_samples','baseline_sample_length','baselineRange'}))
                vd = updateBaseline(vd);
            elseif any(strcmp(updateVar,{'nStd','consecBins','latencyType','latencyRange',...
                    'responsiveAlpha','cusumBaseline','cusumDelta','cusum_nStd',...
                    'min_responsive_trials','responsive_nStd_over_baseline',...
                    'min_fraction_responsive_trials','min_responsive_length'}))
                vd = updateLatency(vd);
            end
            
        end
        function vd = checkClass(vd,responsiveOnly,varargin)
            if isempty(varargin)
                vd.misclassified = nan(1,vd.nCells);
                c = 1;
            else
                c = varargin{1};
            end
            
            if responsiveOnly
                cell_ks = vd.responsiveCells(find(vd.responsiveCells>=c,1):end);
            else
                usableCells = find(vd.usable);
                cell_ks = usableCells(find(usableCells>=c,1):end);
            end
            figure('units','normalized','outerposition',[0 0 1 1]);
            hold on
            for c = cell_ks
                subplot(2,1,1)
                vd.plotRaster(c)
                subplot(2,1,2)
                vd.plotFR(c);
                title([vd.batNum{c} '-' strjoin(strsplit(vd.cellInfo{c},'_'),'-')]);
                commandwindow
                class = input('misclassified?');
                if class == 1 || class == 0
                    vd.misclassified(c) = class;
                else
                    vd.misclassified(c) = Inf;
                    break
                end
                clf
            end
            
            
        end
        function plotFR(vd,cell_k,resp_trial_flag)
            tRange = [min(vd.time) max(vd.time)];
            if strcmp(vd.call_echo,'call')
                lineColor = 'r';
            elseif strcmp(vd.call_echo,'echo')
                lineColor = 'b';
            end
            
            if nargin < 3
                resp_trial_flag = false;
            end
            
            if resp_trial_flag == 1
                resp_trials = get_responsive_trials(vd,cell_k,vd.responsive_nStd_over_baseline,vd.latencyRange);
                
                used_calls = find(vd.usedCalls{cell_k});
                if any(resp_trials)
                    fr = median(vd.trialFR{cell_k}(used_calls(resp_trials),:));
                    frDev = mad(vd.trialFR{cell_k}(used_calls(resp_trials),:))/sqrt(sum(resp_trials));
                else
                   fr = zeros(1,length(vd.time));
                   frDev = zeros(1,length(vd.time));
                end
            elseif resp_trial_flag == 2
                used_calls = vd.usedCalls{cell_k};
                fr = mean(vd.trialFR{cell_k}(used_calls,:));
                frDev = std(vd.trialFR{cell_k}(used_calls,:))/sqrt(sum(used_calls));
            else
                fr = vd.avgFR{cell_k};
                frDev = vd.devFR{cell_k}';
            end
            
            boundedline(vd.time,fr,frDev,lineColor);
            if vd.avgBaseline(cell_k) - vd.devBaseline(cell_k) > 0
                plot(tRange,repmat(vd.avgBaseline(cell_k) - vd.devBaseline(cell_k),1,2),'k.-','LineWidth',2)
                plot(tRange,repmat(vd.avgBaseline(cell_k) - 2*vd.devBaseline(cell_k),1,2),'k-','LineWidth',1)
            else
                plot(tRange,zeros(1,2),'k.-','LineWidth',2)
            end
            plot(tRange,repmat(vd.avgBaseline(cell_k) + vd.devBaseline(cell_k),1,2),'k.-','LineWidth',2)
            plot(tRange,repmat(vd.avgBaseline(cell_k) + 2*vd.devBaseline(cell_k),1,2),'k-','LineWidth',1)
            plot(repmat(vd.latency(cell_k),1,2),get(gca,'ylim'),'k:','LineWidth',2);
            xlim(tRange);
            
            xlabel('Time (s)');
            ylabel('Firing Rate (Hz)');
            commandwindow;
        end
        function plotRaster(vd,cell_k,varargin)
            hold on
            tRange = [min(vd.time) max(vd.time)];
            responsive_order = false;
            if ~isempty(varargin)
                if numel(varargin{1}) > 1
                    order = varargin{1};
                    if size(order,1) ~= 1
                        order = num2cell(order');
                    else
                        order = num2cell(order);
                    end
                else
                    resp_trials = get_responsive_trials(vd,cell_k,vd.responsive_nStd_over_baseline,vd.latencyRange);
                    used_calls = find(vd.usedCalls{cell_k});
                    order = 1:length(used_calls);
                    order(resp_trials) = 1:sum(resp_trials);
                    order(~resp_trials) = sum(resp_trials)+1:length(used_calls);
                    order = num2cell(order);
                    responsive_order = true;
                end
                
            else
                order = num2cell(1:sum(vd.usedCalls{cell_k}));
            end
            cellfun(@(trial,r) plot(r,repmat(trial,length(r)),'k.','MarkerSize',12),order,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}))
            xlim(tRange);
            ylim([0 sum(vd.usedCalls{cell_k})+1]);
            xlabel('Time (s)');
            ylabel('Trial #');
            if responsive_order
               plot(get(gca,'XLim'),repmat(sum(resp_trials)+1.5,1,2),'g--') 
            end
        end
        function plotCalls(vd,cData,cell_k,varargin)
            
            hold on
            [callIdx, callTrain] = callDataIdx(vd,cData,cell_k);
            
            if ~isempty(varargin)
                order = varargin{1};
                callIdx = callIdx(order);
                callTrain = callTrain(order);
            end
            
            
            if strcmp(vd.onset_or_offset,'onset')
                hCalls = plot([zeros(length(callIdx),1) cData.callLength(callIdx)]',repmat(1:length(callIdx),2,1),'r','LineWidth',7);
            elseif strcmp(vd.onset_or_offset,'offset')
                hCalls = plot([-cData.callLength(callIdx) zeros(length(callIdx),1)]',repmat(1:length(callIdx),2,1),'r','LineWidth',7);
            end
            for k = 1:length(callIdx)
                hCalls(k).Color(4) = 0.33;
                hTrain = plot(callTrain(k).relative_callPos',repmat(k,2,size(callTrain(k).relative_callPos,1)),'g','LineWidth',7);
                for l = 1:length(hTrain)
                    hTrain(l).Color(4) = 0.33;
                end
            end
            
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
        function [idx, callTrain] = callDataIdx(vd,cData,cell_k)
            
            idx = find(ismember(cData.callID, vd.callNum{cell_k}(vd.usedCalls{cell_k})));     
            exp_day_call_idx = (cData.expDay == vd.expDay(cell_k)) & strcmp(cData.batNum,vd.batNum{cell_k});

            if nargout > 1
                callTrain_idx = cell(1,sum(vd.usedCalls{cell_k}));
                callTrain_callPos = cell(1,sum(vd.usedCalls{cell_k}));
                for k = 1:length(idx)
                    if strcmp(vd.onset_or_offset,'onset')
                        callTrain_idx{k} = find(cData.callPos(:,1) - cData.callPos(idx(k),1) > 0 &...
                            cData.callPos(:,1) - cData.callPos(idx(k),1) < vd.spikeRange(2) & ...
                            exp_day_call_idx);
                        callTrain_callPos{k} = vertcat(cData.callPos(callTrain_idx{k},:) - cData.callPos(idx(k),1));
                    elseif strcmp(vd.onset_or_offset,'offset')
                        callTrain_idx{k} = find(cData.callPos(:,1) - cData.callPos(idx(k),1) < 0 &...
                            cData.callPos(:,1) - cData.callPos(idx(k),1) > vd.spikeRange(1) & ...
                            exp_day_call_idx);
                        callTrain_callPos{k} = vertcat(cData.callPos(callTrain_idx{k},:) - cData.callPos(idx(k),2));
                    end
                end
                
                callTrain = struct('relative_callPos',callTrain_callPos,'idx',callTrain_idx);
                
            end
        end
        function boutWF = getBoutWF(vd,cData,cell_k,call_k)
            [idx,callTrain] = callDataIdx(vd,cData,cell_k);
            bout_call_idx = [callTrain(call_k).idx];
            bout_call_pos = [0 cData.callLength(idx(call_k)); callTrain(call_k).relative_callPos];
            boutWF = cData.callWF{idx(call_k)}';
            for k = 1:length(bout_call_idx)
                boutWF = [boutWF zeros(1,round(cData.fs*(bout_call_pos(k+1,1)-bout_call_pos(k,2)))) cData.callWF{bout_call_idx(k)}'];
            end
        end
        function [boutCallWF,bout_t,csc,csc_t,callSpikes] = get_call_bout_csc_spikes(vd,cData,bcData,cell_ks,call_k,csc_nums)
            call_offset = 0.2;
            filter_cutoff_frequencies=[600 6000];
            b = find(strcmp(vd.batNums,vd.batNum{cell_ks(1)}));
            expDate = vd.cellInfo{cell_ks(1)}(1:strfind(vd.cellInfo{cell_ks(1)},vd.tetrodeStr)-1);
            base_dir = [vd.baseDirs{b} 'bat' vd.batNum{cell_ks(1)} filesep 'neurologger_recording' expDate filesep];
            trial_used_calls = find(vd.usedCalls{cell_ks(1)});
            [~,cut_call_data,audio2nlg] = get_cell_info(vd,cell_ks(1));
            [~,callTrain] = callDataIdx(vd,cData,cell_ks(1));
            cp = repmat(cut_call_data(trial_used_calls(call_k)).corrected_callpos(1),1,2) + 1e3.*[-call_offset max(callTrain(call_k).relative_callPos(:))+call_offset];
            
            callSpikes = cell(1,length(cell_ks));
            for k = 1:length(cell_ks)
                callSpikes{k} = vd.callSpikes{cell_ks(k)}{trial_used_calls(call_k)};
            end
            
            csc_t = cell(1,length(csc_nums));
            csc = cell(1,length(csc_nums));
            filter_loaded = false;
            for k = 1:length(csc_nums)
                tsData = load([base_dir 'nlxformat' filesep 'CSC' num2str(csc_nums(k)) '.mat']);
                if ~filter_loaded
                    sampling_freq=1/(tsData.sampling_period_usec/1e6);
                    [b,a]=butter(6,filter_cutoff_frequencies/(sampling_freq/2),'bandpass');
                    filter_loaded = true;
                end
                nSamp = length(tsData.AD_count_int16);
                timestamps_usec = get_timestamps_for_Nlg_voltage_all_samples(nSamp,tsData.indices_of_first_samples,tsData.timestamps_of_first_samples_usec,tsData.sampling_period_usec);
                timestamps_msec = 1e-3*timestamps_usec - audio2nlg.first_nlg_pulse_time;
                [~,csc_idx] = inRange(timestamps_msec,cp);
                csc{k} = filtfilt(b,a,double(tsData.AD_count_to_uV_factor*tsData.AD_count_int16(csc_idx)));
                csc_t{k} = linspace(-call_offset,1e-3*diff(cp)-call_offset,length(csc{k}));
            end
            
            boutIdx = find(cellfun(@(x) any(ismember(cData.callID(callTrain(call_k).idx),x)),bcData.boutCalls));
            boutCallWF = get_bout_WF(bcData,cData,boutIdx);
            boutCallWF = [zeros(1,call_offset*cData.fs) boutCallWF zeros(1,call_offset*cData.fs)];
            bout_t = linspace(-call_offset,-call_offset+length(boutCallWF)/250e3,length(boutCallWF));
            
        end
        function [spikeCorr,p,spikes1,spikes2] = scCorr(vd,cell_1,cell_2,period)
            if (strcmp(vd.batNum{cell_1},vd.batNum{cell_2}) && vd.expDay(cell_1) == vd.expDay(cell_2))
                [~,idx1,idx2] = intersect(vd.callNum{cell_1},vd.callNum{cell_2});
                callSpikes1 = vd.callSpikes{cell_1}(idx1);
                callSpikes2 = vd.callSpikes{cell_2}(idx2);
                switch period
                    case 'prePost'
                        spikes1 = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,vd.postCall]))/(vd.preCall+vd.postCall),callSpikes1);
                        spikes2 = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,vd.postCall]))/(vd.preCall+vd.postCall),callSpikes2);
                        
                    case 'pre'
                        spikes1 = cellfun(@(spikes) length(inRange(spikes,[0,-vd.preCall]))/vd.preCall,callSpikes1);
                        spikes2 = cellfun(@(spikes) length(inRange(spikes,[0,-vd.preCall]))/vd.preCall,callSpikes2);
                        
                    case 'post'
                        spikes1 = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,callSpikes1);
                        spikes2 = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,callSpikes2);
                        
                    case 'during'
                        spikes1 = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,callSpikes1,num2cell(vd.callLength{cell_1}(idx2)));
                        spikes2 = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,callSpikes2,num2cell(vd.callLength{cell_2}(idx2)));
                        
                    case 'base'
                        spikes1 = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/abs(diff(vd.baselineRange)),callSpikes1);
                        spikes2 = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/abs(diff(vd.baselineRange)),callSpikes2);
                        
                        
                end
                [spikeCorr, p] = corrcoef(spikes1,spikes2);
                spikeCorr = spikeCorr(1,2);
                p = p(1,2);
                
            else
                disp('No spike count correlation for neurons on different sessions');
                return
            end
        end
        function [leadingCalls, boutCalls, nonBoutCalls] = FR_by_bout_calls(vd,cData,thresh,cell_ks,timeWin)
            
            if nargin < 4
                cell_ks = vd.responsiveCells;
                timeWin = [-vd.preCall vd.postCall];
            elseif nargin < 5
                timeWin = [-vd.preCall vd.postCall];
            end
            
            callFR = cell(3,vd.nCells);
            
            for cell_k = cell_ks
                idx = ismember(cData.callID, vd.callNum{cell_k});
                callPos = cData.callPos(idx,:);
                pre_interCallinterval = [Inf; (callPos(2:end,1) - callPos(1:end-1,2))]; % list of how far current call is from last call
                post_interCallinterval = abs([(callPos(1:end-1,2) - callPos(2:end,1)); Inf]); % list of how far current call is from last call
                
                leading_call_idx = pre_interCallinterval>thresh & post_interCallinterval<thresh;
                inBout_idx = ~leading_call_idx & (pre_interCallinterval<thresh | post_interCallinterval<thresh);
                outBout_idx = post_interCallinterval>=thresh & pre_interCallinterval>=thresh;
                
                k = 1;
                for call_idx = {leading_call_idx,inBout_idx,outBout_idx}
                    idx = call_idx{1};
                    if sum(idx) >= vd.minCalls 
                        trialSpikes = vd.callSpikes{cell_k};
                        all_call_length = num2cell(vd.callLength{cell_k});
                        assert(length(trialSpikes) == (sum(leading_call_idx) + sum(inBout_idx) + sum(outBout_idx)));
                        callFR{k,cell_k} = cellfun(@(spikes,cLength) length(inRange(spikes,timeWin + [0 cLength]))/range((timeWin + [0 cLength])),trialSpikes(idx),all_call_length(idx));
                    end
                    k = k + 1;
                end
            end
            leadingCalls = callFR(1,:);
            boutCalls = callFR(2,:);
            nonBoutCalls = callFR(3,:);
            
        end
        function [low_FR, high_FR] = FR_by_audio_feature(vd,cData,varName,thresh,cell_ks,timeWin)
            
             if nargin < 5
                cell_ks = vd.responsiveCells;
                timeWin = [-vd.preCall vd.postCall];
            elseif nargin < 6
                timeWin = [-vd.preCall vd.postCall];
            end
            n_used_cells = length(cell_ks);
            
            low_FR = cell(1,n_used_cells);
            high_FR = cell(1,n_used_cells);
            
            min_used_calls = 10;
            
            for k = 1:n_used_cells
                cell_k = cell_ks(k);
                
                nTrial = length(vd.callSpikes{cell_k});
                idx = ismember(cData.callID, vd.callNum{cell_k});
                low_idx = find(cData.(varName)(idx)<thresh);
                high_idx = find(cData.(varName)(idx)>=thresh);
                if length(low_idx) >= min_used_calls && length(high_idx) >= min_used_calls
                    trialSpikes = vd.callSpikes{cell_k};
                    cLength = vd.callLength{cell_k};
                    fr = zeros(1,nTrial);
                    for tt = 1:nTrial
                        fr(tt) = length(inRange(trialSpikes{tt},timeWin + [0 cLength(tt)]))/range(timeWin + [0 cLength(tt)]);
                    end
                    low_FR{k} = fr(low_idx);
                    high_FR{k} = fr(high_idx);
                end
            end
        end

        function cells = simulResp(vd,varargin)
            if ~isempty(varargin)
                resp_cells_by_bat = varargin{1};
            else
                resp_cells_by_bat = vd.responsive_cells_by_bat;
            end
            k = 1;
            cells = {};
            for b = 1:3
                respCells = resp_cells_by_bat;
                for cell_k = resp_cells_by_bat
                    respCells = setdiff(respCells,cell_k);
                    idx = find(vd.expDay(respCells) == vd.expDay(cell_k));
                    if ~isempty(idx)
                        cells{k} = [cell_k respCells(idx)];
                        respCells = setdiff(respCells,respCells(idx));
                        k = k+1;
                    end
                end
            end
            
        end
        function [callFR, out_of_call_FR] = getCallFR(vd,cData,cell_k,callOffset,bootFlag)
            [~, callTrain] = callDataIdx(vd,cData,cell_k);
            trialIdx = vd.usedCalls{cell_k};
            if nargin == 3 
                callOffset = zeros(1,2);
                bootFlag = true;
            elseif nargin == 4
                bootFlag = true;
            end
            
            out_of_call_FR = zeros(1,length(callTrain));
            callFR = zeros(1,length(callTrain));
            used_call_spikes = vd.callSpikes{cell_k}(trialIdx);
            used_call_lengths = vd.callLength{cell_k}(trialIdx);
            for bout_k = 1:length(callTrain)
                callPos = [0 used_call_lengths(bout_k); callTrain(bout_k).relative_callPos];
                callPos = callPos + repmat(callOffset,size(callPos,1),1);
                spikes = used_call_spikes{bout_k};
                spike_counts = histcounts(spikes,[vd.time(1)-vd.dT vd.time]);
                call_idx = false(1,length(vd.time));
                for call_k = 1:size(callPos,1)
                    [~,idx] = inRange(vd.time,callPos(call_k,:));
                    call_idx = call_idx | idx;
                end
                call_time = sum(call_idx)*vd.dT;
                callFR(bout_k) = sum(spike_counts(call_idx))/call_time;
                
                if bootFlag
                    bootFR = nan(1,vd.nBoot);
                    bout_call_lengths = diff(callPos,[],2);
                    maxT = min(vd.spikeRange(end),max(callPos(:,2)));
                    minT = 0;
                    for b = 1:vd.nBoot
                        boot_call_idx = false(1,length(vd.time));
                        for boot_call_k = 1:size(callPos,1)
                            call_start = (maxT-minT)*rand+minT;
                            boot_call_pos = [call_start call_start + bout_call_lengths(boot_call_k)];
                            [~,idx] = inRange(vd.time,boot_call_pos);
                            boot_call_idx = boot_call_idx | idx;
                        end
                        boot_call_time = sum(boot_call_idx)*vd.dT;
                        bootFR(b) = sum(spike_counts(boot_call_idx))/boot_call_time;
                    end
                    out_of_call_FR(bout_k) = mean(bootFR);
                else
                    spikeBounds = [-1 min(vd.spikeRange(end),max(callPos(:,2))+1)];
                    [~,spike_bounds_idx] = inRange(vd.time,spikeBounds);
                    non_call_idx = ~call_idx & spike_bounds_idx;
                    non_call_time = sum(non_call_idx)*vd.dT;                    
                    out_of_call_FR(bout_k) = sum(spike_counts(non_call_idx))/non_call_time;
                end
            end
        end
        function [resp_trials, resp_callIDs, non_resp_callIDs] = get_responsive_trials(vd,cell_k,n_std_over_baseline,respRange,trialIdx)

            if nargin < 3
                n_std_over_baseline = 2;
                respRange = [-1 1];
                trialIdx = vd.usedCalls{cell_k};
            elseif nargin < 4
                respRange = [-1 1];
                trialIdx = vd.usedCalls{cell_k};
            elseif nargin < 5
                trialIdx = vd.usedCalls{cell_k};
            end
            
            [~,time_idx] = inRange(vd.time,respRange);
            fr = vd.trialFR{cell_k}(trialIdx,:);
            resp_trials = any(fr(:,time_idx) > vd.avgBaseline(cell_k) + n_std_over_baseline*vd.devBaseline(cell_k),2);
            callIDs = vd.callNum{cell_k}(trialIdx);
            resp_callIDs = callIDs(resp_trials);
            non_resp_callIDs = callIDs(~resp_trials);
        end
        function call_info = get_call_bhv_info(vd,cell_k)
            if strcmp(vd.batNum{cell_k},vd.batNums{4})
                
                expDate = vd.cellInfo{cell_k}(1:strfind(vd.cellInfo{cell_k},vd.tetrodeStr)-1);
                audioDir = [vd.baseDirs{4} 'bat' vd.batNum{cell_k} filesep 'neurologger_recording' expDate '\audio\ch1\'];
                
                try
                    s = load([audioDir 'juv_call_info_' vd.batNum{cell_k} '_' vd.call_echo '.mat']);
                    call_info = s.call_info;
                catch err
                    disp(err)
                    keyboard
                end
                
            else
                bhv_file_dir = 'E:\ephys\juvenile_recording\bhvFiles\';
                exp_datestr = datestr(vd.expDay(cell_k),'yyyymmdd');
                bhvFName = [bhv_file_dir strjoin({'juv_call_info',vd.batNum{cell_k},exp_datestr},'_')];
                
                try
                    s = load(bhvFName);
                    try
                        call_info = s.juv_call_info;
                    catch
                        try
                            call_info = s.call_info;
                        catch err
                            disp(err)
                            keyboard
                        end
                    end
                catch
                    disp('couldn''t find bhv file')
                    call_info = [];
                end
            end
                
        end
        function bhvs = callID2AudioFile(vd,cData,cell_k)
            
            call_info = get_call_bhv_info(vd,cell_k);
            
            if isempty(call_info)
               bhvs = [];
               return
            end
            
            if strcmp(vd.batNum{cell_k},vd.batNums{4})
                
                bhvs = call_info;
                
            else
                
                callID = vd.callNum{cell_k};
                
                juv_call_times = arrayfun(@(x) datetime(x.AudioFile(strfind(x.AudioFile,'T')+1:strfind(x.AudioFile,'_')-1),'InputFormat','yyMMddHHmmss'),call_info);
                bhvs = struct('AudioFile',[],'juvCall',[],'echoCall',[],'VideoFile',[],'behaviors',[]);
                
                for ID_k = 1:length(callID)
                    fName = cData('callID',callID(ID_k)).fName;
                    [fPath,fName] = fileparts(fName{:});
                    fPath = strsplit(fPath,filesep);
                    fPath = [strjoin(fPath(1:end-1),filesep) filesep];
                    file_call_pos = cData('callID',callID(ID_k)).file_call_pos;
                    callTime = datetime(fName(strfind(fName,'T')+1:strfind(fName,'_')-1),'InputFormat','yyMMddHHmmss') + seconds(file_call_pos(1)/cData.fs);
                    [~,idx] = min(abs(juv_call_times-callTime));
                    
                    callData = audioread([fPath 'ch2' filesep call_info(idx).AudioFile]);
                    callWF = cData('callID',callID(ID_k)).callWF;
                    match = strfind(callData',callWF');
                    if ~isempty(match)  %#ok<STREMP>
                        bhvs(ID_k) = call_info(idx);
                    else
                        bhvs(ID_k) = struct('AudioFile',[],'juvCall',[],'echoCall',[],'VideoFile',[],'behaviors',[]);
                    end
                end
                
                callID = num2cell(callID);
                [bhvs.callID] = deal(callID{:});
            end
        end
        function bhvIdx = find_call_bhv(vd,cData,cell_k,varargin)

            bhv = callID2AudioFile(vd,cData,cell_k);
            nBhvSeq = length(bhv);
            bhvIdx = true(1,nBhvSeq);
            
            for s = 1:nBhvSeq
                if ~isempty(bhv(s).behaviors)
                    nBhv = sum(~cellfun(@isempty,bhv(s).behaviors));
                    incBhv = true(1,nBhv);
                    for b = 1:nBhv
                        behaviorString = bhv(s).behaviors{b};
                        behaviorStringSplit = strsplit(behaviorString,'-');
                        for subs_k = 1:2:length(varargin)
                            switch varargin{subs_k}
                                
                                case 'juvAdult'
                                    incBhv(b) = incBhv(b) && strcmpi(behaviorStringSplit{1},varargin{subs_k+1});
                                    
                                case 'contact'
                                    incBhv(b) = incBhv(b) && strcmpi(behaviorStringSplit{2},varargin{subs_k+1});
                                    
                                case 'direction'
                                    incBhv(b) = incBhv(b) && strcmpi(behaviorStringSplit{3},varargin{subs_k+1});
                                    
                                case 'bhv'
                                    incBhv(b) = incBhv(b) && any(strcmpi(behaviorStringSplit{4},varargin{subs_k+1}));
                            end
                        end
                    end
                    bhvIdx(s) = any(incBhv);
                    
                else
                    bhvIdx(s) = false;
                end
            end
            
        end
        function [bhv_trial_idx, non_bhv_trial_idx] = get_calls_by_bhv(vd,cData,cell_k,selectedBhvs)
            if nargin < 4
                selectedBhvs = {'bite','strike','claw'}; % select aggresive behaviors
            end
            
            try
                bhvIdx = find_call_bhv(vd,cData,cell_k,'bhv',selectedBhvs,'juvAdult','Juvenile');
            catch
                bhv_trial_idx = [];
                non_bhv_trial_idx = [];
                return
            end
            
            bhv_trial_idx = vd.usedCalls{cell_k} & bhvIdx;
            non_bhv_trial_idx = vd.usedCalls{cell_k} & ~bhvIdx;
        end
        function [p,bhv_fr,non_bhv_fr] = FR_by_behavior(vd,cData,cell_k,minCalls,selectedBhvs,responseRange)
            
            if nargin < 4
                minCalls = vd.minCalls;
                selectedBhvs = {'bite','strike','claw'}; % select aggresive behaviors
                responseRange = [-0.1 0.1];
            elseif nargin < 5
                selectedBhvs = {'bite','strike','claw'}; % select aggresive behaviors
                responseRange = [-0.1 0.1];
            elseif nargin < 6
                responseRange = [-0.1 0.1];
            end
            try
                [bhv_trial_idx, non_bhv_trial_idx] = get_calls_by_bhv(vd,cData,cell_k,selectedBhvs);
            catch
                p = NaN;
                bhv_fr = NaN;
                non_bhv_fr = NaN;
                return
            end
            
            [~,t_idx] = inRange(vd.time,responseRange);
            fr = vd.trialFR{cell_k};
            bhv_fr = mean(fr(bhv_trial_idx,t_idx),2);
            non_bhv_fr = mean(fr(non_bhv_trial_idx,t_idx),2);
            if length(bhv_fr) >= minCalls && length(non_bhv_fr) >= minCalls
                p = ranksum(bhv_fr, non_bhv_fr);
            else
                p = NaN;
            end
        end
        function latency = trialLatency(vd,cell_k,trialIdx)
            if nargin < 3
               trialIdx = vd.usedCalls{cell_k}; 
            end
            switch vd.latencyType
                case 'cusum'
                    [~,~,~,latency]  = cusumResponsive(vd,cell_k,trialIdx);
                case 'trialBased'
                    [~, ~, ~, latency]  = trialBasedResponsive(vd,cell_k,trialIdx);
                case 'std_above_baseline'
                    [~, ~, ~, latency]  = std_above_baseline_responsive(vd,cell_k,trialIdx);
                    
            end
        end
        function total_call_length = get_total_call_length(vd,cData,cell_k,call_offset)
            [idx, callTrain] = callDataIdx(vd,cData,cell_k);
            used_call_length = cData.callLength(idx);
            
            total_call_length = 0;
            for k = 1:length(callTrain)
                if isempty(callTrain(k).relative_callPos)
                    total_call_length = total_call_length + used_call_length(k) + call_offset;
                else
                    total_call_length = total_call_length + max(callTrain(k).relative_callPos(:)) + call_offset;
                end
            end
        end
        function [total_session_length, sessionBounds] = get_total_session_length(vd,cell_k)
            
            expDate = vd.cellInfo{cell_k}(1:strfind(vd.cellInfo{cell_k},vd.tetrodeStr)-1);
            b = strcmp(vd.batNums,vd.batNum{cell_k});
            nlx_dir = [vd.baseDirs{b} 'bat' vd.batNum{cell_k} filesep 'neurologger_recording' expDate '\nlxformat\'];
            events = load([nlx_dir 'EVENTS.mat']);
            [stabilityBounds, ~, audio2nlg] = getCellInfo(vd,b,cell_k);
            stabilityBounds = 1e-3*stabilityBounds;
            
            session_strings = {'start_communication','end_communication'};
            alternate_session_strings = {'end_playback','Stopped recording'};
            sessionTime = zeros(1,length(session_strings));
            sessionBounds = zeros(1,length(session_strings));
            for s = 1:length(session_strings)
                sessionIdx = cellfun(@(x) contains(x,session_strings{s}),events.event_types_and_details);
                if ~sum(sessionIdx)
                    sessionIdx = cellfun(@(x) contains(x,alternate_session_strings{s}),events.event_types_and_details);
                end
                if s == 2 && ~sum(sessionIdx)
                    sessionIdx = length(events.event_timestamps_usec)-1;
                end
                sessionTime(s) = 1e-3*(1e-3*events.event_timestamps_usec(sessionIdx) - audio2nlg.first_nlg_pulse_time);
                
                if s == 1
                    
                    sessionBounds(s) = max(sessionTime(s),stabilityBounds(s));
                    
                elseif s == 2
                    
                    sessionBounds(s) = min(sessionTime(s),stabilityBounds(s));
                    
                end
            end
            total_session_length = abs(diff(sessionBounds));
        end
        function [fr_to_session_calls,percent_high_fr_calls,percent_session_calls,cell_ks] = get_fr_to_session_call_ratio(vd,cData,b)
            high_fr_call_info_fnames = dir([vd.baseDirs{b} 'bat' vd.batNums{b} '\**\call_info*high_fr*.mat']);
            percent_session_calls = zeros(1,length(high_fr_call_info_fnames));
            percent_high_fr_calls = zeros(1,length(high_fr_call_info_fnames));
            cell_ks = zeros(1,length(high_fr_call_info_fnames));
            for k = 1:length(high_fr_call_info_fnames)
                idx = strfind(high_fr_call_info_fnames(k).name,'_20');
                cell_info_str = high_fr_call_info_fnames(k).name(idx+1:end-4);
                cell_ks(k) = find(strcmp(vd.cellInfo,cell_info_str));
                total_call_length = get_total_call_length(vd,cData,cell_ks(k),2);
                total_session_length = get_total_session_length(vd,cell_ks(k));
                percent_session_calls(k) = 100*total_call_length/total_session_length;
                s = load(fullfile(high_fr_call_info_fnames(k).folder,high_fr_call_info_fnames(k).name));
                call_info = s.call_info;
                voc_events = arrayfun(@(y) any(cellfun(@(x) ~isempty(x) && contains(x,'Vocalization'),y.behaviors)),call_info);
                percent_high_fr_calls(k) = 100*sum(voc_events)/length(voc_events);
            end
            
            fr_to_session_calls = percent_high_fr_calls./percent_session_calls;
            fr_to_session_calls(fr_to_session_calls==0) = min(fr_to_session_calls);
        end
        function [stabilityBounds, cut_call_data, audio2nlg, ttDir, success] = get_cell_info(vd,cell_k)
            b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
            [stabilityBounds, cut_call_data, audio2nlg, ttDir, success] = getCellInfo(vd,b,cell_k);
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
        function vd = update_usedCalls_by_callData(vd,cData,selecFun)
            
            for cell_k = find(vd.usable)
                try
                    
                    b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
                    [stabilityBounds, cut_call_data] = getCellInfo(vd,b,cell_k);
                    callNums = vd.callNum{cell_k};
                    selectedCalls = selecFun(cData,callNums)';
                    usedCallNums = get_used_calls(stabilityBounds,cut_call_data,1e3*vd.spikeRange,vd.onset_or_offset,selectedCalls);
                    
                    usedCallIdx = usedCallNums;
                    vd.usedCalls{cell_k} = usedCallIdx;
                catch err
                    keyboard
                end
                
            end
            
            vd = updateUsedCells(vd);
            
        end
        function export_spike_data(vd,output_dir)
            
            lastProgress = 0;
            usedCells = false(1,vd.nCells);
            for cell_k = 1:vd.nCells
                b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
                
                switch vd.sortingMetric
                    case 'sortingQuality'
                        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
                    case 'ratio_and_distance'
                        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
                end
                if wellSorted
                    timestamps = getSpikes(vd,b,cell_k);
                    fName = [output_dir vd.batNum{cell_k} '_' vd.cellInfo{cell_k} '.csv'];
                    dlmwrite(fName,timestamps,'delimiter',',','precision','%3f');
                    usedCells(cell_k) = true;
                end
                
                progress = 100*(cell_k/vd.nCells);
                
                if mod(progress,10) < mod(lastProgress,10)
                    fprintf('%d %% of cells processed\n',round(progress));
                end
                
                lastProgress = progress;
                
                
            end
            
        end
        function timestamps = getSpikes(vd,b,cell_k)
            
            [stabilityBounds, ~, audio2nlg, ttDir] = getCellInfo(vd,b,cell_k);
            
            switch vd.expType
                
                case 'adult'
                    
                    timestamps = csvread([vd.spike_data_dir vd.batNum{cell_k} '_' vd.cellInfo{cell_k} '.csv']);
                    if size(timestamps,1) ~= 1
                        timestamps = timestamps';
                    end
                    
                case 'juvenile'
                    timestamps = Nlx2MatSpike(ttDir,[1 0 0 0 0],0,1,[]); % load sorted cell data
                    timestamps = 1e-3*timestamps - audio2nlg.first_nlg_pulse_time; % convert to ms and align to first TTL pulse on the NLG
                    timestamps = inRange(timestamps,stabilityBounds);
            end
        end
        function export_all_call_spike_data(vd,output_dir)
            lastProgress = 0;
            usedCells = false(1,vd.nCells);
            for cell_k = 1:vd.nCells
                b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
                
                switch vd.sortingMetric
                    case 'sortingQuality'
                        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
                    case 'ratio_and_distance'
                        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
                end
                if wellSorted
                    timestamps = getSpikes(vd,b,cell_k);
                    [stabilityBounds, cut_call_data] = getCellInfo(vd,b,cell_k);
                    fName = [output_dir vd.batNum{cell_k} '_' vd.cellInfo{cell_k} '.mat'];
                    all_call_info = struct('timestamps',timestamps,'cut_call_data',cut_call_data,'stabilityBounds',stabilityBounds);
                    save(fName,'-v7.3','-struct','all_call_info')
                    usedCells(cell_k) = true;
                end
                
                progress = 100*(cell_k/vd.nCells);
                
                if mod(progress,10) < mod(lastProgress,10)
                    fprintf('%d %% of cells processed\n',round(progress));
                end
                
                lastProgress = progress;
                
                
            end
        end
    end
    
end

function [stabilityBounds, cut_call_data, audio2nlg, ttDir, success] = getCellInfo(vd,b,cell_k)

cellInfo = vd.cellInfo{cell_k};
baseDir = vd.baseDirs{b};
batNum = vd.batNums{b};



switch vd.expType
    
    case 'adult'
        
        expDate = cellInfo(1:strfind(cellInfo,vd.tetrodeStr)-1);
        audioDir = [baseDir 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
        stabilityBounds = [-Inf Inf];
        
        
        try
            
            s = load([audioDir 'cut_call_data.mat']);
            cut_call_data = s.cut_call_data;
            cut_call_data = cut_call_data(~[cut_call_data.noise]);
            if isempty(cut_call_data)
                audio2nlg = [];
                ttDir = [];
                success = true;
                return
            end
            batIdx = unique(cellfun(@(call) find(cellfun(@(bNum) strcmp(bNum,vd.batNums{b}),call)),{cut_call_data.batNum}));
            
            if length(batIdx) == 1
                callpos = horzcat(cut_call_data.corrected_callpos);
                callpos = callpos(batIdx,:);
                [cut_call_data.corrected_callpos] = deal(callpos{:});
            else
                keyboard
            end
            call_info_fname = dir([audioDir 'call_info_*_' vd.call_echo '_' expDate '.mat']);
            if length(call_info_fname) > 1
                keyboard
            end
            s = load(fullfile(call_info_fname.folder,call_info_fname.name));
            call_info = s.call_info;
            
            assert(all([cut_call_data.uniqueID] == [call_info.callID]));
            
            bat_calls = cellfun(@(x) ischar(x{1}) && contains(x,batNum),{call_info.behaviors});
            cut_call_data = cut_call_data(bat_calls);
            
            audio2nlg = [];
            ttDir = [];
            
        catch err
            
            disp(err)
            keyboard
            success = false; %input('continue?');
            
            if ~success
                
                audio2nlg = [];
                ttDir = [];
                return
            end
            
        end
        success = true;
        
    case 'juvenile'
        
        cluster_num_idx = strfind(cellInfo,vd.clusterStr)+length(vd.clusterStr);
        tt_num_idx = strfind(cellInfo,vd.tetrodeStr )+length(vd.tetrodeStr );
        
        expDate = cellInfo(1:strfind(cellInfo,vd.tetrodeStr)-1);
        tetrodeNum = str2double(cellInfo(tt_num_idx));
        cellNum = str2double(cellInfo(cluster_num_idx:end));
        call_echo = vd.call_echo;
        
        if cellNum < 10
            sorted_cell_string = [vd.clusterStr '0'];
        else
            sorted_cell_string = vd.clusterStr;
        end
        
        tetrode = [expDate vd.tetrodeStr num2str(tetrodeNum) sorted_cell_string num2str(cellNum)]; % build tetrode/cell string
        ttDir = [vd.spike_data_dir 'bat' batNum filesep expDate filesep tetrode '.ntt']; % filename for cell of interest
        audioDir = [baseDir 'bat' batNum filesep 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
        s = load([vd.analysisDir 'cell_stability_info']);
        cell_stability_info = s.cell_stability_info;
        idx = find(strcmp(tetrode,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
        stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
        
        
        try
            
            if strcmp(call_echo,'call') % are we looking at calls or echolocation clicks?
                s = load([audioDir 'cut_call_data.mat']);
                cut_call_data = s.cut_call_data;
                cut_call_data = cut_call_data(~[cut_call_data.noise]);
                
            elseif strcmp(call_echo,'echo')
                s = load([audioDir 'cut_echo_data.mat']);
                cut_call_data = s.cut_call_data;
                cut_call_data = cut_call_data(~[cut_call_data.noise]);
                
                s = load([audioDir 'juv_call_info_' vd.batNums{b} '_echo.mat']);
                call_info = s.call_info;
                
                assert(all([cut_call_data.uniqueID] == [call_info.callID]));
                
                echo_calls = strcmp({call_info.echoCall},'juvEcho');
                cut_call_data = cut_call_data(echo_calls);
            end
            
        catch err
            
            disp(err)
            keyboard
            
            success = input('continue?');
            
            if ~success
                return
            end
            
        end
        
        audio2nlg = load([audioDir 'audio2nlg_fit.mat']); % load fit data to sync audio to nlg data
        stabilityBounds = stabilityBounds - audio2nlg.first_nlg_pulse_time;
        success = true;
end

end

function [callSpikes, callNum, callLength] = get_call_spikes(timestamps,stabilityBounds,cut_call_data,spikeRange,call_onset_offset,timeWarp,warp_call_length)

nCalls = length(cut_call_data); % total # of calls, excluding any errors by cutting procedure

call_k = 1; % counter for all used calls

callSpikes = cell(1,nCalls); % initialize cell to contain spike times relative to call onset
callNum = nan(1,nCalls);
callLength = nan(1,nCalls);

for call = 1:nCalls % iterate through all the calls within this .wav file
    cp = cut_call_data(call).corrected_callpos; % time of onset and offset of the call of interest on this iteration
    callposOffset = [cp(1)+spikeRange(1), cp(2)+spikeRange(2)]; % total time window around call of interest
    
    if ~any(isnan(cp)) % check that the TTL alignment worked for this call
        if (callposOffset(1) > stabilityBounds(1) && callposOffset(2) < stabilityBounds(2)) % check that this call is within timeframe that this cell is stable
            
            relativeCP = cp(strcmp({'onset','offset'},call_onset_offset)); % choose beginning or end of call
            
            timestampsCall = inRange(timestamps,callposOffset); % all spikes occuring within that time window
            callSpikes{call_k} = 1e-3*(timestampsCall - relativeCP); % align spikes to call onset/offset, convert to sec, and store
            callNum(call_k) = cut_call_data(call).uniqueID;
            callLength(call_k) = 1e-3*diff(cp);
            if timeWarp
                callSpikes{call_k} = callSpikes{call_k} * warp_call_length/callLength(call_k);
            end
            call_k = call_k + 1; % increment total calls
        end
        
    end
end

%now truncate the list of calls/clicks to how many calls/clicks actually occurred
nCalls = call_k - 1; % determine the actual number of calls
callSpikes = callSpikes(1:nCalls); % truncate list of spikes relative to call onset to length of n_calls
callNum = callNum(1:nCalls);
callLength = callLength(1:nCalls);
end

function usedCalls = get_used_calls(stabilityBounds,cut_call_data,spikeRange,call_onset_offset,varargin)

nCalls = length(cut_call_data); % total # of calls

if ~isempty(varargin)
    selectedCalls = varargin{1};
else
    selectedCalls = [cut_call_data.uniqueID];
end

all_call_idx = ismember([cut_call_data.uniqueID],selectedCalls);
all_call_times = [cut_call_data(all_call_idx).corrected_callpos];

usedCalls = true(1,nCalls);
call_k = 1;
rec_k = 1;
for call = 1:nCalls % iterate through all the calls within this .wav file
    cp = cut_call_data(call).corrected_callpos; % time of onset and offset of the call of interest on this iteration
    callposOffset = [cp(1)+spikeRange(1), cp(2)+spikeRange(2)]; % total time window around call of interest
    
    if ~any(isnan(cp)) % check that the TTL alignment worked for this call
        if (callposOffset(1) > stabilityBounds(1) && callposOffset(2) < stabilityBounds(2)) % check that this call is within timeframe that this cell is stable
            if (strcmp(call_onset_offset,'onset') && any(all_call_times > (cp(1) + spikeRange(1)) & all_call_times < cp(1))) || ... skip any calls within a call train
               (strcmp(call_onset_offset,'offset') && any(all_call_times < (cp(2) + spikeRange(2)) & all_call_times > cp(2))) || ...
               ~(ismember(cut_call_data(call).uniqueID,selectedCalls))    
                usedCalls(call_k) = false;
            end
            call_k = call_k + 1;
            rec_k = rec_k + 1;
            
        else
            rec_k = rec_k + 1;
        end
    end
end

nCalls = call_k - 1; % determine the actual number of calls
usedCalls = usedCalls(1:nCalls);

end

function [baselineMean, baselineStd] = calculate_baseline(vd,cell_k,timestamps,cut_call_data)

switch vd.baselineMethod
    case 'preCall'
        baselineLength = abs(diff(vd.baselineRange));
        nTrial = sum(vd.usedCalls{cell_k});
        allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        usedSpikes = inRange(allSpikes,vd.baselineRange);
        baselineMean = length(usedSpikes)/nTrial/baselineLength;
%         baselineMean = median(cellfun(@(x) length(inRange(x,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k})));
        baselineStd = std(cellfun(@(x) length(inRange(x,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k})));
    case 'randomSamples'
        if nargin < 3
            b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
            [~, cut_call_data, ~, ~] = getCellInfo(vd,b,cell_k);
            timestamps = getSpikes(vd,b,cell_k);
        end
        
        if strcmp(vd.expType,'adult')
            timestamps = timestamps - timestamps(1);
        end
        
        baselineSamples = cell(1,vd.n_baseline_samples);
        baselineLength = vd.baseline_sample_length*1e3;
        all_call_times = [cut_call_data.corrected_callpos];
        n_sample = 0;
        max_t = max(timestamps) - baselineLength;
        while n_sample < vd.n_baseline_samples
            samp_t = max_t*rand;
            if ~any(all_call_times > (samp_t - baselineLength/2) & all_call_times < (samp_t + baselineLength/2))
                n_sample = n_sample + 1;
                baselineSamples{n_sample} = inRange(timestamps,[(samp_t - baselineLength/2) (samp_t + baselineLength/2)]);
            end
        end
        baselineMean = length([baselineSamples{:}])/vd.baseline_sample_length/vd.n_baseline_samples;
        baselineStd = std(cellfun(@(x) length(x)/vd.baseline_sample_length,baselineSamples));
        
    case 'sessionFR'
        
        if nargin < 3
            b = find(strcmp(vd.batNums,vd.batNum{cell_k}));
            [~, cut_call_data, ~, ~] = getCellInfo(vd,b,cell_k);
            timestamps = getSpikes(vd,b,cell_k);
        end
        
        if strcmp(vd.expType,'adult')
            timestamps = timestamps - timestamps(1);
        end
        
        [fr, t] = sskernel(1e-6*timestamps(timestamps>0),linspace(0,max(1e-6*timestamps),1e5));
        fr = fr/1e3;
        t = t*1e6;
        call_t = true(1,length(t));
        baselineLength = vd.baseline_sample_length*1e3;
        
        for k = 1:length(cut_call_data)
            call_t(t > (cut_call_data(k).corrected_callpos(1) - baselineLength) & t < (cut_call_data(k).corrected_callpos(2) + baselineLength)) = false;
        end
        
        baselineMean = mean(fr);
        baselineStd = std(fr);
end

end

function usable = checkUsability(vd,cell_k)

% usable = all(cellfun(@(x) length(x(x > vd.callRange(1) & x < vd.callRange(2))),vd.callSpikes{cell_k})>=vd.minSpikes) && length(vd.callSpikes{cell_k}) >= vd.minCalls;
allSpikes = vd.callSpikes{cell_k}(vd.usedCalls{cell_k});
allSpikes = [allSpikes{:}];
prePostSpikes = inRange(allSpikes,[-vd.preCall,vd.postCall]);
nTrial = sum(vd.usedCalls{cell_k});
switch vd.sortingMetric
    case 'sortingQuality'
        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
    case 'ratio_and_distance'
        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
end
usable =  (length(prePostSpikes)/sum(vd.usedCalls{cell_k})) >= vd.minSpikes &&...
    nTrial >= vd.minCalls &&...
    wellSorted;

end

function [spike_rate_mean, spike_rate_dev, spikeRate, spikeTrain, bw] = frKernelEstimate(vd,cell_k)

time = vd.callRange(1):vd.dT:vd.callRange(2);

nTrial = length(vd.callSpikes{cell_k});
nT = length(time);
spikeRate = zeros(nTrial,nT);
spikeTrain = zeros(nTrial,nT);

switch vd.kernelType
    case 'manual'
        
        if isnan(vd.constantBW)
            bw = kde([vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}]);
        else
            bw = vd.constantBW;
        end
        
        for tt = 1:nTrial
            kSpike = 1;
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[time(1),time(end)]);
            if ~isempty(usedSpikes)
                n_used_spikes = length(usedSpikes);
                spike_rate_trial = zeros(n_used_spikes,nT);
                spikeTrain(tt,:) = histcounts(usedSpikes,[time time(end)+vd.dT]);
                for s = 1:n_used_spikes
                    kernel = Gauss(time,usedSpikes(s),bw.^2);
                    spike_rate_trial(kSpike,:) = kernel*vd.dT;
                    kSpike = kSpike + 1;
                end
                spikeRate(tt,:) = sum(spike_rate_trial,1)/vd.dT;
            end
            
        end
        
        
    case 'fixedOpt'
        bw = nan(nTrial,1);
        for tt = 1:nTrial
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[time(1),time(end)]);
            if ~isempty(usedSpikes)
                spikeTrain(tt,:) = histcounts(usedSpikes,[time time(end)+vd.dT]);
                [spikeRate(tt,:), ~, bw(tt)]= sskernel(usedSpikes,time,vd.frBandwidthRange);
            end
        end
        
    case 'varOpt'
        bw = nan(nTrial,nT);
        for tt = 1:nTrial
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[time(1),time(end)]);
            if ~isempty(usedSpikes)
                spikeTrain(tt,:) = histcounts(usedSpikes,[time time(end)+vd.dT]);
                [spikeRate(tt,:), ~, bw(tt,:)]= ssvkernel(usedSpikes,time,vd.frBandwidthRange);
            end
        end
        
       
end

[spike_rate_mean, spike_rate_dev] = calculateAvgFR(vd,cell_k,spikeRate);

[~,idx] = inRange(time,[time(1)+vd.smoothingOffset,time(end)-vd.smoothingOffset]);
spikeRate = spikeRate(:,idx);
spikeTrain = spikeTrain(:,idx);

end

function [spike_rate_mean, spike_rate_dev] = calculateAvgFR(vd,cell_k,spikeRate)

time = vd.callRange(1):vd.dT:vd.callRange(2);

switch vd.kernelType
    case 'manual'
        
        spike_rate_mean = mean(spikeRate(vd.usedCalls{cell_k},:));
        spike_rate_dev = std(spikeRate(vd.usedCalls{cell_k},:))/sqrt(size(spikeRate,1))';
        spike_rate_dev = repmat(spike_rate_dev,2,1);
        
    case 'fixedOpt'
        
        allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        [spike_rate_mean,~,~,~,~,spike_rate_dev] = sskernel(allSpikes,time);
        spike_rate_mean = spike_rate_mean/sum(vd.usedCalls{cell_k});
        
    case 'varOpt'
        
        allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        [spike_rate_mean,~,~,~,~,spike_rate_dev] = ssvkernel(allSpikes,time);
        spike_rate_mean = spike_rate_mean/sum(vd.usedCalls{cell_k});
        
end

[~,idx] = inRange(time,[time(1)+vd.smoothingOffset,time(end)-vd.smoothingOffset]);
spike_rate_mean = spike_rate_mean(idx);
spike_rate_dev = spike_rate_dev(:,idx);

end

function y = Gauss(t,tS,w)

y = (1./sqrt(2*pi*w)) * exp( -(((t - tS).^2)./(2*w)));

end

function [latency, respType, respValency, respStrength] = calculateLatency(vd,cell_k)


if any(strcmp(vd.latencyType,{'pre','post','during','prePost','peri'}))
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
    
elseif strcmp(vd.latencyType,'cusum')
    [respType, respValency, respStrength, latency] = cusumResponsive(vd,cell_k);
    if isempty(inRange(latency,vd.latencyRange))
        latency = NaN;
        respType = NaN;
        respValency = NaN;
        respStrength = NaN;
    end
    
elseif strcmp(vd.latencyType,'trialBased')
    [respType, respValency, respStrength, latency]  = trialBasedResponsive(vd,cell_k);
elseif strcmp(vd.latencyType,'std_above_baseline')
    [respType, respValency, respStrength, latency]  = std_above_baseline_responsive(vd,cell_k);

end


end

function [respType, respValency, respStrength, latency]  = cusumResponsive(vd,cell_k,trialIdx)

usedCalls = find(vd.usedCalls{cell_k});
if nargin < 3
    trialIdx = 1:sum(vd.usedCalls{cell_k});
end


spikeTrain = sum(vd.trial_spike_train{cell_k}(usedCalls(trialIdx),:),1);
idx = find(vd.time>=vd.cusumStart);
muHat = mean(spikeTrain(vd.time<vd.cusumBaseline));
SInc = zeros(1,length(idx)+1);
SDec = zeros(1,length(idx)+1);

mStd_Inc = vd.cusumDelta;
mStd_Dec = 1/vd.cusumDelta;
nStd = vd.cusum_nStd;

for tt = 1:length(idx)
    SInc(tt+1) = max(0,SInc(tt) + spikeTrain(idx(tt))*log(mStd_Inc) + (1-mStd_Inc)*muHat);
    SDec(tt+1) = max(0,SDec(tt) + spikeTrain(idx(tt))*log(mStd_Dec) + (1-mStd_Dec)*muHat);
end

if any([sum(SDec>nStd)>1, sum(SInc>nStd)>1])
    
    
    if sum(SInc>nStd)>1 && ~(sum(SDec>nStd)>1)
        
        respValency = 1;
        latencyIdx = idx(find(SInc>nStd,1,'first'));
        
    elseif sum(SDec>nStd)>1 && ~(sum(SInc>nStd)>1 )
        
        respValency = -1;
        latencyIdx = idx(find(SDec>nStd,1,'first'));
        
        
    elseif sum(SDec>nStd)>1 && sum(SInc>nStd)>1
        
        latencyIdx_Inc = idx(find(SInc>nStd,1,'first'));
        latencyIdx_Dec = idx(find(SDec>nStd,1,'first'));
        
        respStrength_Inc = abs(vd.avgFR{cell_k}(latencyIdx_Inc) - vd.avgBaseline(cell_k));
        respStrength_Dec = abs(vd.avgFR{cell_k}(latencyIdx_Dec) - vd.avgBaseline(cell_k));
        
        if respStrength_Inc > respStrength_Dec
            respValency = 1;
            latencyIdx = latencyIdx_Inc;
        else
            respValency = -1;
            latencyIdx = latencyIdx_Dec;
        end
        
    end
    
    latency = vd.time(latencyIdx);
    if latency < 0
        respType = 'pre';
    else
        respType = 'post';
    end
    
    respStrength = abs(vd.avgFR{cell_k}(latencyIdx) - vd.avgBaseline(cell_k));
    
else
    latency = NaN;
    respType = NaN;
    respValency = NaN;
    respStrength = NaN;
end

end

function [responsive, respType, respValency, respStrength] = periCallResponsive(vd,cell_k)

baselineLength = diff(vd.baselineRange);
baseSpikes = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));

switch vd.latencyType
    
    case 'during'
        callSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}),num2cell(vd.callLength{cell_k}(vd.usedCalls{cell_k})));
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'during';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'pre'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'pre';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'post'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'post';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'prePost'
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));
        postCallSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));
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
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}));
        postCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[callLength,callLength + vd.postCall]))/vd.postCall,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}),num2cell(vd.callLength{cell_k}(vd.usedCalls{cell_k})));
        duringCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k}),num2cell(vd.callLength{cell_k}(vd.usedCalls{cell_k})));
        responsivePre = ttest(baseSpikes,preCallSpikes,'Alpha',vd.responsiveAlpha);
        responsivePost = ttest(baseSpikes,postCallSpikes,'Alpha',vd.responsiveAlpha);
        responsiveDuring = ttest(baseSpikes,duringCallSpikes,'Alpha',vd.responsiveAlpha);
        if any(isnan([responsivePre , responsivePost , responsiveDuring]))
            responsive = false;
        else
            responsive = responsivePre | responsivePost | responsiveDuring;
        end
            
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

function [respType, respValency, respStrength, latency]  = trialBasedResponsive(vd,cell_k,trialIdx)

if nargin < 3
    trialIdx = vd.usedCalls{cell_k};
end

[respType, respValency, respStrength, latency] = deal(NaN);
%%
resp_trials = get_responsive_trials(vd,cell_k,vd.responsive_nStd_over_baseline,vd.latencyRange,trialIdx);

if sum(resp_trials) >= vd.min_responsive_trials && sum(resp_trials)/length(resp_trials) > vd.min_fraction_responsive_trials

    [latency_time,time_idx] = inRange(vd.time,vd.latencyRange);
    
    used_calls = find(trialIdx);
    
    fr = vd.trialFR{cell_k}(used_calls(resp_trials),:);
    latency_time_fr = fr(:,time_idx);
    
    baseline_avg = max(vd.avgBaseline(cell_k),1);
    
    over_thresh = median(latency_time_fr,1) - median(mad(fr,1,1))/sqrt(size(fr,1)) > baseline_avg + vd.responsive_nStd_over_baseline*vd.devBaseline(cell_k);
    
    if any(over_thresh)
        thresh_crossing = find(diff(over_thresh) == 1) + 1;
        crossing_length = zeros(1,length(thresh_crossing));
        
        k = 1;
        for tc = thresh_crossing
            if ~isempty(find(diff(over_thresh(tc:end)) == -1,1))
                crossing_length(k) = find(diff(over_thresh(tc:end)) == -1,1);
            else
                crossing_length(k) = length(time_idx) - tc;
            end
            k = k + 1;
        end
        
        wins = [thresh_crossing' (thresh_crossing  + crossing_length)'];
        merged_wins = merge_wins(wins,1/vd.dT,vd.min_responsive_length);
        thresh_crossing = merged_wins(:,1)';
        crossing_length = diff(merged_wins,[],2)';
        
        if any(crossing_length*vd.dT >= vd.min_responsive_length)
            thresh_crossing_idx = thresh_crossing(find(crossing_length*vd.dT >= vd.min_responsive_length,1,'first'));
            latency = latency_time(thresh_crossing_idx);
            if latency < 0
                respType = 'pre';
            else
                respType = 'post';
            end
            respValency = 1;
            respStrength = median(latency_time_fr(over_thresh)) - vd.avgBaseline(cell_k);
        end
        
    end
    %%
end

end

function [respType, respValency, respStrength, latency]  = std_above_baseline_responsive(vd,cell_k)

if nargin < 3
    trialIdx = vd.usedCalls{cell_k};
end

[respType, respValency, respStrength, latency] = deal(NaN);
[latency_time,time_idx] = inRange(vd.time,vd.latencyRange);

fr = mean(vd.trialFR{cell_k}(trialIdx,:));
frDev = std(vd.trialFR{cell_k}(trialIdx,:))/sqrt(sum(trialIdx));

thresh = vd.avgBaseline(cell_k) + vd.responsive_nStd_over_baseline*vd.devBaseline(cell_k);

latency_time_fr = fr(:,time_idx);

over_thresh = latency_time_fr - frDev(time_idx) > thresh;

if any(over_thresh)
    thresh_crossing = find(diff(over_thresh) == 1) + 1;
    crossing_length = zeros(1,length(thresh_crossing));
    
    k = 1;
    for tc = thresh_crossing
        if ~isempty(find(diff(over_thresh(tc:end)) == -1,1))
            crossing_length(k) = find(diff(over_thresh(tc:end)) == -1,1);
        else
            crossing_length(k) = length(time_idx) - tc;
        end
        k = k + 1;
    end
    
    wins = [thresh_crossing' (thresh_crossing  + crossing_length)'];
    merged_wins = merge_wins(wins,1/vd.dT,vd.min_responsive_length);
    thresh_crossing = merged_wins(:,1)';
    crossing_length = diff(merged_wins,[],2)';
    
    if any(crossing_length*vd.dT >= vd.min_responsive_length)
        thresh_crossing_idx = thresh_crossing(find(crossing_length*vd.dT >= vd.min_responsive_length,1,'first'));
        latency = latency_time(thresh_crossing_idx);
        if latency < 0
            respType = 'pre';
        else
            respType = 'post';
        end
        respValency = 1;
        respStrength = median(latency_time_fr(over_thresh)) - vd.avgBaseline(cell_k);
    end
    
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
            [vd.avgFR{cell_k}, vd.devFR{cell_k}, vd.trialFR{cell_k},...
                vd.trial_spike_train{cell_k}, vd.frBandwidth{cell_k}] =  frKernelEstimate(vd,cell_k);
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

function vd = updateUsedCells(vd)

vd = updateUsability(vd);

for cell_k = find(vd.usable)
    [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k,vd.trialFR{cell_k});

end
vd = updateBaseline(vd);


end
