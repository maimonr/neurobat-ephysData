classdef vocalData < ephysData
    properties(SetAccess = public)
        spikeRange = [-1 1]
        callRange = [-1 1]
        baselineRange = [-1 0]
        frBandwidthRange = 2e-3:5e-3:1e-1
        kernelType = 'manual'
        dT = 5e-3
        minCalls = 15
        minSpikes = 1e-3
        nStd = 2
        consecBins = 15;
        nBoot = 1000
        sortingMetric = 'sortingQuality'
        sortingThreshold = 2
        latencyType = 'cusum'
        latencyRange = [-1 1]
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
        cusumBaseline = -1.5
        cusumStart = -1.5
        waveformSamples = 32;
        timeWarp = false
        warp_call_length = 0.1
        baselineMethod = 'preCall'
        n_baseline_samples = 1e2
        baseline_sample_length = 5
        exclude_neighboring_calls = true
        operant_reward_status = 'rewardedOnly'
        selectCalls = 'selfCall'
        cellType = 'singleUnit'
        clusterStr = '_SS_'
        tetrodeStr = 'TT'
        onset_or_offset = 'onset'
        
        frBandwidth
        usedCalls
        isolationDistance
        LRatio
        sortingQuality
        daysOld
        batNum
        cellInfo
        expDay
        callSpikes
        callNum
        call_bat_num
        callLength
        trial_spike_train
        stored_trial_fr
        avgFR
        devFR
        avgBaseline
        devBaseline
        latency
        respType
        respValency
        respStrength
        usable
        meanFR
        medianISI
        peak2trough
        spikeWidth
        avg_spike
        duration
        tetrodeNum
        tetrodeDepth
    end
    
    properties (SetAccess = private)
        misclassified
    end
    
    properties (Dependent)
        time
        responsiveCells
        responsiveIdx 
        responsive_cells_by_bat
        nCells
        cellTable
    end
    methods
        function vd = vocalData(varargin)
            
            pnames = {'callType','onset_or_offset','expType','vd_update','selectCells','selectCalls','cellType','minCalls','operant_reward_status'};
            dflts  = {'call','onset','juvenile',[],[],'selfCall','singleUnit',15,'rewardedOnly'};
            [callType,onset_or_offset,expType,vd_update,selectCells,selectCalls,cellType,minCalls,operant_reward_status] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            vd = vd@ephysData(expType);
            
            vd.callType = callType;
            vd.onset_or_offset = onset_or_offset;
            vd.selectCalls = selectCalls;
            vd.cellType = cellType;
            vd.minCalls = minCalls;
            vd.operant_reward_status = operant_reward_status;
            
            if ~isempty(vd_update)
                vd = vd_update;
            end
            
            use_select_cells = ~isempty(selectCells);
            
            vd.expDay = datetime([],[],[]);
            
            cell_k = 1;
            lastProgress = 0;
            
            if any(strcmp(vd.expType{1},{'adult','adult_operant'}))
                
                nBats = length(vd.batNums);
                spike_file_names = dir([vd.spike_data_dir{1} '*.csv']);
                single_unit_idx = arrayfun(@(fName) contains(fName.name,vd.clusterStr),spike_file_names);
                
                switch vd.cellType
                    case 'singleUnit'
                        spike_file_names = spike_file_names(single_unit_idx);
                        
                        sortingInfo = load([vd.analysisDir{1} 'sortingInfo.mat']);
                        sortingInfo = sortingInfo.sortingInfo;
                        cellInfo_regexp_str = '\d{8}_TT\d_SS_\d{2}';
                        
                    case 'multiUnit'
                        spike_file_names = spike_file_names(~single_unit_idx);
                        sortingInfo = [];
                        cellInfo_regexp_str = '\d{8}_TT\d';
                end
                
                nCells_to_process = length(spike_file_names);
                
                [vd.isolationDistance, vd.LRatio, vd.sortingQuality, vd.avgBaseline,...
                    vd.devBaseline, vd.respValency, vd.respStrength, vd.meanFR,...
                    vd.medianISI, vd.peak2trough, vd.spikeWidth, vd.duration,...
                    vd.daysOld, vd.latency, vd.respType, vd.tetrodeNum, vd.tetrodeDepth]  = deal(nan(1,nCells_to_process));
                [vd.batNum, vd.cellInfo, vd.callSpikes, vd.callNum, vd.frBandwidth,...
                    vd.callLength, vd.stored_trial_fr, vd.avgFR, vd.devFR,...
                    vd.trial_spike_train, vd.usedCalls, vd.call_bat_num] = deal(cell(1,nCells_to_process));
                vd.respType = num2cell(vd.respType);
                vd.usable = false(1,nCells_to_process);
                vd.avg_spike = nan(nCells_to_process,vd.waveformSamples);
                
                tt_depth_format = '%{MM/dd/yy}D %f %f %f %f';
                ttString = 'TT';
                
                for b = 1:nBats
                    bat_spike_file_names = spike_file_names(arrayfun(@(x) contains(x.name,vd.batNums{b}),spike_file_names));
                    
                    try
                        tetrodeDepths = readtable(fullfile(vd.baseDirs{b},'tetrode_depths',['tetrode_depths_' vd.batNums{b} '.csv']),'format',tt_depth_format);
                        use_tetrode_depths = true;
                    catch
                        use_tetrode_depths = false;
                    end
                    
                    for d = 1:length(bat_spike_file_names)
                        
                        cellInfo = regexp(bat_spike_file_names(d).name,cellInfo_regexp_str,'match');
                        cellInfo = cellInfo{1};
                        if use_select_cells && ~any(strcmp(cellInfo,selectCells))
                            cell_k = cell_k + 1;
                            continue
                        end
                        
                        vd.batNum{cell_k} = vd.batNums{b};
                        vd.cellInfo{cell_k} = cellInfo;
                        vd.tetrodeNum(cell_k) = str2double(cellInfo(strfind(cellInfo,ttString)+length(ttString)));
                        
                        % get sorting quality for this cell
                        if strcmp(vd.cellType,'singleUnit') 
                            sortingInfo_idx = strcmp({sortingInfo.cellInfo},cellInfo) & strcmp({sortingInfo.batNum},vd.batNum{cell_k});
                            if strcmp(vd.expType{1},'adult')
                                vd.isolationDistance(cell_k) = sortingInfo(sortingInfo_idx).isolationDistance;
                                vd.LRatio(cell_k) = sortingInfo(sortingInfo_idx).LRatio;
                                vd.sortingQuality(cell_k) = sortingInfo(sortingInfo_idx).sortingQuality;
                            elseif strcmp(vd.expType{1},'adult_operant')
                                
                                if strcmp(vd.callType,'call') % sorting quality for this experiment is calculated for each session
                                    sorting_info_str = 'communication';
                                elseif any(strcmp(vd.callType,{'operant','operant_reward'}))
                                    sorting_info_str = 'operant';
                                end
                                vd.isolationDistance(cell_k) = sortingInfo(sortingInfo_idx).isolationDistance.(sorting_info_str);
                                vd.LRatio(cell_k) = sortingInfo(sortingInfo_idx).LRatio.(sorting_info_str);
                                vd.sortingQuality(cell_k) = sortingInfo(sortingInfo_idx).sortingQuality.(sorting_info_str);
                            end
                        end
                        
                        vd.expDay(cell_k) = datetime(cellInfo(1:8),'InputFormat',vd.dateFormat);
                        if use_tetrode_depths
                            try
                                vd.tetrodeDepth(cell_k) = tetrodeDepths{tetrodeDepths.Date==vd.expDay(cell_k),vd.tetrodeNum(cell_k)+1};
                            catch
                                vd.tetrodeDepth(cell_k) = NaN;
                            end
                        else
                            vd.tetrodeDepth(cell_k) = NaN;
                        end
                        
                        [stabilityBounds, cut_call_data] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
                        
                        if isempty(cut_call_data)
                            cell_k = cell_k + 1;
                            continue
                        end
                        
                        timestamps = getSpikes(vd,cell_k);
                        
                        [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.call_bat_num{cell_k}, vd.callLength{cell_k},vd.usedCalls{cell_k}] = ...
                            get_used_call_spikes(vd,stabilityBounds,cut_call_data,'timestamps',timestamps);
                        
                        vd.daysOld(cell_k) = days(vd.expDay(cell_k) - vd.birthDates{b});
                        
                        if checkUsability(vd,cell_k)
                            vd.usable(cell_k) = true;
                            [vd.stored_trial_fr{cell_k},vd.trial_spike_train{cell_k},...
                                vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
                            [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k);
                            [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k,timestamps,cut_call_data);
                            [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                        end
                        
                        progress = 100*(cell_k/nCells_to_process);
                        
                        if mod(progress,10) < mod(lastProgress,10)
                            fprintf('%d %% of cells processed\n',round(progress));
                        end
                        
                        lastProgress = progress;
                        
                        cell_k = cell_k + 1;
                    end
                end
                
            elseif strcmp(vd.expType{1},'adult_wujie')
                
                spike_file_names = dir([vd.spike_data_dir{1} '*.csv']);
                spike_file_names = spike_file_names(arrayfun(@(x) contains(x.name,vd.batNums),spike_file_names));
                spike_file_name_dlm = '_';
                spike_file_name_format = 'bbbbb_yyyymmddTTt_SS_ss';
                nBats = length(vd.batNums);
                nCells_to_process = length(spike_file_names);
                
                [vd.avgBaseline, vd.devBaseline, vd.respValency,...
                    vd.respStrength, vd.meanFR, vd.latency, vd.respType,...
                    vd.tetrodeNum, vd.tetrodeDepth]  = deal(nan(1,vd.nCells));
                [vd.batNum, vd.cellInfo, vd.callSpikes, vd.callNum, vd.frBandwidth,...
                    vd.callLength, vd.stored_trial_fr, vd.avgFR, vd.devFR,...
                    vd.trial_spike_train, vd.usedCalls] = deal(cell(1,vd.nCells));
                vd.respType = num2cell(vd.respType);
                vd.usable = false(1,vd.nCells);
                for b = 1:nBats
                    spike_file_names = dir([vd.spike_data_dir{1} '*.csv']);
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
                        [stabilityBounds, cut_call_data, ~, ~, ~, success] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
                        timestamps = csvread([vd.spike_data_dir{1} spike_file_names(d).name]);
                        
                        if size(timestamps,1) ~= 1
                            timestamps = timestamps';
                        end
                        
                        [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.call_bat_num{cell_k}, vd.callLength{cell_k},vd.usedCalls{cell_k}] = ...
                            get_used_call_spikes(vd,stabilityBounds,cut_call_data,'timestamps',timestamps);
                        
                        vd.sortingQuality(cell_k) = vd.sortingThreshold;
                        if checkUsability(vd,cell_k) && success
                            vd.usable(cell_k) = true;
                            [vd.stored_trial_fr{cell_k},vd.trial_spike_train{cell_k},...
                                vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
                            [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k);
                            [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k,timestamps,cut_call_data);
                            [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                        end
                        
                        progress = 100*(cell_k/nCells_to_process);
                        
                        if mod(progress,10) < mod(lastProgress,10)
                            fprintf('%d %% of cells processed\n',round(progress));
                        end
                        
                        lastProgress = progress;
                        
                        cell_k = cell_k + 1;
                    end
                end
                
                
            elseif strcmp(vd.expType{1},'juvenile')
                addpath('C:\Users\phyllo\Documents\Maimon\ephys\scripts\experimentation_scripts\Wujie\MatlabImportExport_v6.0.0\')
                nBats = length(vd.batNums);
                sortedCells = load([vd.analysisDir{1} 'sortedCells.mat']);
                sortingInfo = load([vd.analysisDir{1} 'sortingInfo.mat']);
                sortingInfo = sortingInfo.sortingInfo ;
                nCells_to_process = length(sortedCells.batNumList);
                
                [vd.isolationDistance, vd.LRatio, vd.sortingQuality, vd.avgBaseline,...
                    vd.devBaseline, vd.respValency, vd.respStrength, vd.meanFR,...
                    vd.medianISI, vd.peak2trough, vd.spikeWidth, vd.duration,...
                    vd.daysOld, vd.latency, vd.respType, vd.tetrodeNum, vd.tetrodeDepth]  = deal(nan(1,nCells_to_process));
                [vd.batNum, vd.cellInfo, vd.callSpikes, vd.callNum, vd.frBandwidth,...
                    vd.callLength, vd.stored_trial_fr, vd.avgFR, vd.devFR,...
                    vd.trial_spike_train, vd.usedCalls] = deal(cell(1,nCells_to_process));
                vd.respType = num2cell(vd.respType);
                vd.usable = false(1,nCells_to_process);
                vd.avg_spike = nan(nCells_to_process,vd.waveformSamples);
                
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
                        
                        [stabilityBounds, cut_call_data, ~, ttDir] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
                        timestamps = getSpikes(vd,cell_k);
                        
                        [vd.callSpikes{cell_k}, vd.callNum{cell_k}, vd.call_bat_num{cell_k}, vd.callLength{cell_k}, vd.usedCalls{cell_k}] = ...
                            get_used_call_spikes(vd,stabilityBounds,cut_call_data,'timestamps',timestamps);
                        
                        [vd.meanFR(cell_k), vd.medianISI(cell_k), vd.peak2trough(cell_k),...
                            vd.spikeWidth(cell_k), vd.avg_spike(cell_k,:), vd.duration(cell_k)] = get_spike_stats(timestamps,ttDir);
                        
                        vd.daysOld(cell_k) = days(vd.expDay(cell_k) - vd.birthDates{b});
                        
                        if checkUsability(vd,cell_k)
                            vd.usable(cell_k) = true;
                            [vd.stored_trial_fr{cell_k},vd.trial_spike_train{cell_k},...
                                vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
                            [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k);
                            [vd.avgBaseline(cell_k), vd.devBaseline(cell_k)] = calculate_baseline(vd,cell_k,timestamps,cut_call_data);
                            [vd.latency(cell_k), vd.respType{cell_k}, vd.respValency(cell_k), vd.respStrength(cell_k)] = calculateLatency(vd,cell_k);
                        end
                        
                        progress = 100*(cell_k/nCells_to_process);
                        
                        if mod(progress,10) < mod(lastProgress,10)
                            fprintf('%d %% of cells processed\n',round(progress));
                        end
                        
                        lastProgress = progress;
                        
                        cell_k = cell_k + 1;
                    end
                end
            end
        end
        function n = numArgumentsFromSubscript(~,~,~)
            n = 1;
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
                                    case 'sortingQuality'
                                        cellIdx = cellIdx & vd.sortingQuality<=S(1).subs{idx+1};
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
                                    case 'cellTable'
                                        cellIdx = cellIdx & ismember(vd.cellTable,S(1).subs{idx+1})';
                                    otherwise
                                        disp('indexing variable not recognized')
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
                                        disp(err)
                                        return
                                    end
                                end
                                
                            else
                                varargout = {vd.(S(2).subs)(cellIdx)};
                            end
                        else
                            try
                                varargout = {builtin('subsref',vd,S)};
                            catch err
                                switch err.message
                                    case 'Too many output arguments.'
                                        builtin('subsref',vd,S);
                                    otherwise
                                        disp('Indexing in VocalData must come in pairs');
                                        return
                                end
                            end
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
                    'responsiveAlpha','cusumBaseline','cusumDelta','cusum_nStd','cusumStart',...
                    'min_responsive_trials','responsive_nStd_over_baseline',...
                    'min_fraction_responsive_trials','min_responsive_length'}))
                vd = updateLatency(vd);
            else
               disp('No associated update function')
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
            
            if ~vd.usable(cell_k)
                return
            end
            
            tRange = [min(vd.time) max(vd.time)];
            if strcmp(vd.callType,'call')
                lineColor = 'r';
            elseif strcmp(vd.callType,'echo')
                lineColor = 'b';
            elseif strcmp(vd.callType,'operant')
                lineColor = 'g';
            elseif strcmp(vd.callType,'operant_reward')
                lineColor = 'cyan';
            end
            
            if nargin < 3
                resp_trial_flag = false;
            end
            
            if resp_trial_flag == 1
                resp_trials = get_responsive_trials(vd,cell_k,vd.responsive_nStd_over_baseline,vd.latencyRange);
                
                used_calls = find(vd.usedCalls{cell_k});
                if any(resp_trials)
                    trialFR = vd.trialFR(cell_k);
                    fr = median(trialFR(used_calls(resp_trials),:));
                    frDev = mad(trialFR(used_calls(resp_trials),:))/sqrt(sum(resp_trials));
                else
                    fr = zeros(1,length(vd.time));
                    frDev = zeros(1,length(vd.time));
                end
            elseif resp_trial_flag == 2
                trialFR = vd.trialFR(cell_k);
                used_calls = vd.usedCalls{cell_k};
                fr = mean(trialFR(used_calls,:));
                frDev = std(trialFR(used_calls,:))/sqrt(sum(used_calls));
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
                else
                    resp_trials = get_responsive_trials(vd,cell_k);
                    used_calls = find(vd.usedCalls{cell_k});
                    order = [used_calls(resp_trials) used_calls(~resp_trials)];
                    responsive_order = true;
                end
                
            else
                order = find(vd.usedCalls{cell_k});
            end
            
            spike_train = vd.trial_spike_train{cell_k};
            spike_train = spike_train(order,:);
            [row,col] = find(spike_train);
            scatter(vd.time(col),row,12,'k','filled')
            
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
        function [fr_sort_idx,fr_time] = plot_fr_heatmap(vd,cell_ks,axisHandle,params)
            
            fr_time = zeros(1,length(cell_ks));
            [t,t_idx] = inRange(vd.time,params.xlims);
            avg_trial_fr = zeros(length(cell_ks),sum(t_idx));
            k = 1;
            for cell_k = cell_ks
                used_calls = vd.usedCalls{cell_k};
                trialFR = vd.trialFR(cell_k);
                avg_trial_fr(k,:) = zscore(nanmean(trialFR(used_calls,t_idx),1));
                if params.resp_valence_flag 
                    avg_trial_fr(k,:) = avg_trial_fr(k,:)*vd.respValency(cell_k);
                end
                [~,idx] = max(avg_trial_fr(k,:));
                fr_time(k) = t(idx);
                k = k + 1;
            end
            
            [fr_time_sort,fr_sort_idx] = sort(fr_time);
            imagesc(axisHandle,t,1:length(cell_ks),avg_trial_fr(fr_sort_idx,:));
            axisHandle.YDir = 'normal';
            caxis(axisHandle,[-1 3])
            axisHandle.XLabel.String = 'Time (s)';
            axisHandle.YLabel.String = 'Cell Number';
            axisHandle.FontSize = params.fontSize;
            colormap(axisHandle,params.cmap)
            yticks_interval = max(50,100*round(length(cell_ks)/3/100));
            axisHandle.YTick = yticks_interval:yticks_interval:length(cell_ks);
            axisHandle.XTick = sort([params.xlims 0]);
            axisHandle.XLim = params.xlims;
            axisHandle.YLim = [1 length(cell_ks)+1];
            axisHandle.Box = 'off';
            axisHandle.TickLength = [0 0];
            axisHandle.PlotBoxAspectRatio = ones(1,3);
            
            if params.plot_fit_line
                lm = fitlm(fr_time_sort,1:length(cell_ks));
                hold(axisHandle,'on')
                x = [min(fr_time) max(fr_time)];
                y = lm.Coefficients.Estimate(1) + lm.Coefficients.Estimate(2)*x;
                plot(x,y,'r--','lineWidth',3)
            end
            
        end
        function [idx, callTrain] = callDataIdx(vd,cData,cell_k)
            
            switch vd.selectCalls
                case 'selfCall'
                    idx = find(ismember(cData.callID, vd.callNum{cell_k}(vd.usedCalls{cell_k})) & strcmp(cData.batNum,vd.batNum(cell_k)) & cData.expDay == vd.expDay(cell_k));
                    exp_day_call_idx = (cData.expDay == vd.expDay(cell_k)) & strcmp(cData.batNum,vd.batNum{cell_k});
                case 'otherCall'
                    idx = find(ismember(cData.callID, vd.callNum{cell_k}(vd.usedCalls{cell_k})) & ~strcmp(cData.batNum,vd.batNum(cell_k)) & cData.expDay == vd.expDay(cell_k));
                    exp_day_call_idx = (cData.expDay == vd.expDay(cell_k)) & ~strcmp(cData.batNum,vd.batNum{cell_k});
                case 'allCall'
                    idx = find(ismember(cData.callID, vd.callNum{cell_k}(vd.usedCalls{cell_k})) & cData.expDay == vd.expDay(cell_k));
                    exp_day_call_idx = (cData.expDay == vd.expDay(cell_k));
                    
            end
            
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
        function [boutCallWF,bout_t,csc,csc_t,callSpikes,batNum,bout_call_nums] = get_call_bout_csc_spikes(vd,cData,bcData,callNum,call_offset,csc_nums)
            cell_ks = find(cellfun(@(cell_call_nums) ismember(callNum,cell_call_nums),vd.callNum));
            batNum = vd.batNum(cell_ks);
            callSpikes = cell(1,length(cell_ks));
            for k = 1:length(cell_ks)
                trial_idx = vd.callNum{cell_ks(k)} == callNum;
                callSpikes{k} = vd.callSpikes{cell_ks(k)}{trial_idx};
            end
            boutSeparation = 1;
            bout_call_nums = get_call_bout_nums(cData,callNum,boutSeparation);
            call_fNames = cData('callID',bout_call_nums).fName;
            file_callPos = cData('callID',bout_call_nums).file_call_pos;
            [~,boutCallWF] = get_bout_WF_manual(bcData,cData,call_fNames,file_callPos([1 end]),call_offset);
  
            bout_call_length = length(boutCallWF)/cData.fs;
            bout_t = linspace(-call_offset,bout_call_length-call_offset,length(boutCallWF));
            
            if nargin > 5
                b = vd.batIdx(cell_ks(1));
                [~,cut_call_data,audio2nlg] = get_cell_info(vd,cell_ks(1));
                [~,callTrain] = callDataIdx(vd,cData,cell_ks(1));
                expDate = vd.cellInfo{cell_ks(1)}(1:strfind(vd.cellInfo{cell_ks(1)},vd.tetrodeStr)-1);
                base_dir = [vd.baseDirs{b} 'bat' vd.batNum{cell_ks(1)} filesep 'neurologger_recording' expDate filesep];
                filter_cutoff_frequencies=[600 6000];
                csc_t = cell(1,length(csc_nums));
                csc = cell(1,length(csc_nums));
                filter_loaded = false;
                call_idx = [cut_call_data.uniqueID] == callNum;
                cp = repmat(cut_call_data(call_idx).corrected_callpos(1),1,2) + 1e3.*[-call_offset max(callTrain(callNum).relative_callPos(:))+call_offset];
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
            else
                csc_t = [];
                csc = [];
            end            
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
        
        function [leadingCalls, boutCalls, nonBoutCalls] = FR_by_bout_calls(vd,cData,thresh,varargin)
            
            pnames = {'cell_ks','timeWin','minCalls'};
            dflts  = {vd.responsiveCells,[-vd.preCall vd.postCall],vd.minCalls};
            [cell_ks,timeWin,minCall] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            if nargin < 4
                cell_ks = vd.responsiveCells;
                timeWin = [-vd.preCall vd.postCall];
            elseif nargin < 5
                timeWin = [-vd.preCall vd.postCall];
            end
            
            callFR = cell(3,length(cell_ks));
            k = 1;
            for cell_k = cell_ks
                idx = ismember(cData.callID, vd.callNum{cell_k}) & cData.expDay == vd.expDay(cell_k);
                callPos = cData.callPos(idx,:);
                pre_interCallinterval = [Inf; (callPos(2:end,1) - callPos(1:end-1,2))]; % list of how far current call is from last call
                post_interCallinterval = abs([(callPos(1:end-1,2) - callPos(2:end,1)); Inf]); % list of how far current call is from last call
                
                leading_call_idx = pre_interCallinterval>thresh & post_interCallinterval<thresh;
                inBout_idx = ~leading_call_idx & (pre_interCallinterval<thresh | post_interCallinterval<thresh);
                outBout_idx = post_interCallinterval>=thresh & pre_interCallinterval>=thresh;
                
                call_type_k = 1;
                for call_idx = {leading_call_idx,inBout_idx,outBout_idx}
                    idx = call_idx{1};
                    if sum(idx) >= minCall
                        trialSpikes = vd.callSpikes{cell_k};
                        all_call_length = num2cell(vd.callLength{cell_k});
                        assert(length(trialSpikes) == (sum(leading_call_idx) + sum(inBout_idx) + sum(outBout_idx)));
                        callFR{call_type_k,k} = cellfun(@(spikes,cLength) length(inRange(spikes,timeWin + [0 cLength]))/range((timeWin + [0 cLength])),trialSpikes(idx),all_call_length(idx));
                    end
                    call_type_k = call_type_k + 1;
                end
                k = k + 1;
            end
            leadingCalls = callFR(1,:);
            boutCalls = callFR(2,:);
            nonBoutCalls = callFR(3,:);
            
        end
        function [low_FR, high_FR, low_idx, high_idx] = FR_by_audio_feature(vd,cData,varName,thresh,cell_ks,timeWin)
            
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
                idx = ismember(cData.callID, vd.callNum{cell_k}) & cData.expDay == vd.expDay(cell_k);
                ICI = [Inf; diff(cData.callPos(:,1))];
                low_idx = find(cData.(varName)(idx)<thresh & ICI(idx) > range(timeWin));
                high_idx = find(cData.(varName)(idx)>=thresh & ICI(idx) > range(timeWin));
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
            used_call_spike_trains = vd.trial_spike_train{cell_k}(trialIdx,:);
            used_call_lengths = vd.callLength{cell_k}(trialIdx);
            for bout_k = 1:length(callTrain)
                callPos = [0 used_call_lengths(bout_k); callTrain(bout_k).relative_callPos];
                callPos = callPos + repmat(callOffset,size(callPos,1),1);
                spike_counts = full(used_call_spike_trains(bout_k,:));
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
        
        function cells = simulResp(vd,varargin)
            if ~isempty(varargin)
                resp_cells_by_bat = varargin{1};
            else
                resp_cells_by_bat = vd.responsive_cells_by_bat;
            end
            k = 1;
            cells = {};
            for b = 1:length(vd.batNums)
                respCells = resp_cells_by_bat{b};
                for cell_k = respCells
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
        function [resp_trials, resp_callIDs, non_resp_callIDs] = get_responsive_trials(vd,cell_k,respRange,trialIdx)
            
            if nargin < 3
                respRange = [-1 1];
                trialIdx = vd.usedCalls{cell_k};
            elseif nargin < 4
                respRange = [-1 1];
                trialIdx = vd.usedCalls{cell_k};
            elseif nargin < 5
                trialIdx = vd.usedCalls{cell_k};
            end
            
            %             [~,time_idx] = inRange(vd.time,respRange);
            %             fr = vd.trialFR{cell_k}(trialIdx,:);
            nSpike = cellfun(@(x) length(inRange(x,respRange))/abs(diff(respRange)),vd.callSpikes{cell_k}(trialIdx));
            p = 1-poisscdf(nSpike,vd.avgBaseline(cell_k));
            resp_trials = p<0.05;
            %             resp_trials = any(fr(:,time_idx) > vd.avgBaseline(cell_k) + n_std_over_baseline*vd.devBaseline(cell_k),2);
            callIDs = vd.callNum{cell_k}(trialIdx);
            resp_callIDs = callIDs(resp_trials);
            non_resp_callIDs = callIDs(~resp_trials);
        end
        function [bat_fr_latency, used_call_idx, cell_ks, manual_selected_cells] = individual_bat_resp(vd,cData,varargin)
            
            pnames = {'cell_ks','respType','minCalls','ICI_limit','plotFlag','timeLims','target_bat_num'};
            dflts  = {vd.responsiveCells,vd.latencyType,vd.minCalls,abs(vd.spikeRange(1)),false,vd.spikeRange,[]};
            [cell_ks,resp_latency_type,minCall,ICI_limit,plotFlag,timeLims,target_bat_num] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            [t,t_idx] = inRange(vd.time,timeLims);
            nCell = length(cell_ks);
            bat_fr_latency = cell(1,nCell);
            used_call_idx = cell(1,nCell);
            manual_selected_cells = nan(1,nCell);
            colors = {'r','g','b','k'};
            k = 1;
            
            for cell_k = cell_ks
                if plotFlag
                    cla
                end
                trialFR = vd.trialFR(cell_k);
                callIDs = vd.callNum{cell_k}';
                call_bat_nums = cData('callID',callIDs,'expDay',vd.expDay(cell_k)).batNum;
                cPos = cData('callID',callIDs,'expDay',vd.expDay(cell_k)).callPos;
                ICI = [Inf; cPos(2:end,1) - cPos(1:end-1,2)];
                select_bat_nums_cat = categorical(setdiff(unique(call_bat_nums(~cellfun(@iscell,call_bat_nums))),[vd.batNum(cell_k) 'unidentified']))';
                idx = ~cellfun(@iscell,call_bat_nums);
                call_bat_nums_cat = cell(length(call_bat_nums),1);
                call_bat_nums_cat(idx) = num2cell(categorical(call_bat_nums(idx)));
                call_bat_nums_cat(~idx) = cellfun(@(x) categorical({x{:}}),call_bat_nums(~idx),'un',0);
                
                if ~isempty(target_bat_num) && ~iscell(target_bat_num)
                    select_bat_nums_cat = {categorical({target_bat_num}),setdiff(select_bat_nums_cat,target_bat_num)};
                elseif ~isempty(target_bat_num) && iscell(target_bat_num)
                    select_bat_nums_cat = {categorical(target_bat_num(k)),setdiff(select_bat_nums_cat,target_bat_num{k})};
                else
                    select_bat_nums_cat = num2cell(select_bat_nums_cat);
                end
                n_select_bats = length(select_bat_nums_cat);
                [bat_fr_latency{k}, used_call_idx{k}, h] = deal(cell(1,n_select_bats));
                for bat_k = 1:n_select_bats 
                    used_call_bat_idx = ICI > ICI_limit & cellfun(@(bNum) any(ismember(bNum,select_bat_nums_cat{bat_k})),call_bat_nums_cat);
                    if sum(used_call_bat_idx)>minCall
                        bat_fr_latency{k}{bat_k} = calculateLatency(vd,cell_k,'latencyType',resp_latency_type,'trialIdx',used_call_bat_idx);
                    end
                    used_call_idx{k}{bat_k} = used_call_bat_idx;
                end
                
                plot_cell_flag = (plotFlag > 0 && plotFlag < 3) || (plotFlag == 3 && any(~isnan([ bat_fr_latency{k}{:}])));
                
                if plot_cell_flag
                    for bat_k = 1:n_select_bats
                        batFR = trialFR(used_call_idx{k}{bat_k},t_idx);
                        if sum(used_call_idx{k}{bat_k})>minCall
                            [~,h{bat_k}] = boundedline(t,mean(batFR),std(batFR)/sqrt(sum(used_call_idx{k}{bat_k})),colors{bat_k},'alpha');
                            h{bat_k} = findobj(h{bat_k});
                        end
                    end
                end
                
                if plot_cell_flag
                    used_idx = cellfun(@sum,used_call_idx{k})>minCall;
                    if any(used_idx)
                        legend_str = cellfun(@(bNum,nCall,lat) cell2mat([strjoin(string(bNum)) '-' num2str(sum(nCall)) '-' num2str(lat)]),select_bat_nums_cat(used_idx),used_call_idx{k}(used_idx),bat_fr_latency{k}(used_idx),'un',0);
                        legend([h{used_idx}],legend_str{:},'Location','northeastoutside')
                        title(sprintf('%s - cell %d',vd.batNum{cell_k},cell_k))
                        axis square
                        xlabel('Time (s)')
                        ylabel('Firing rate (Hz)')
                        set(gca,'FontSize',45);
                        xlim(timeLims)
                        legend box off
                        box off
                        if plotFlag ~= 2
                            selectCell = input('?');
                            if ~isempty(selectCell)
                                manual_selected_cells(k) = selectCell;
                            end
                        end
                    end
                end
                
                k = k + 1;
            end

        end
        
        function call_info = get_call_bhv_info(vd,cell_k)
            if strcmp(vd.batNum{cell_k},vd.batNums{4})
                
                expDate = vd.cellInfo{cell_k}(1:strfind(vd.cellInfo{cell_k},vd.tetrodeStr)-1);
                audioDir = [vd.baseDirs{4} 'bat' vd.batNum{cell_k} filesep 'neurologger_recording' expDate '\audio\ch1\'];
                
                try
                    s = load([audioDir 'juv_call_info_' vd.batNum{cell_k} '_' vd.callType '.mat']);
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
            fr = vd.trialFR(cell_k);
            bhv_fr = mean(fr(bhv_trial_idx,t_idx),2);
            non_bhv_fr = mean(fr(non_bhv_trial_idx,t_idx),2);
            if length(bhv_fr) >= minCalls && length(non_bhv_fr) >= minCalls
                p = ranksum(bhv_fr, non_bhv_fr);
            else
                p = NaN;
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
        function [total_session_length, sessionBounds] = get_total_session_length(vd,cell_k,varargin)
            
            if nargin == 2
                
                [stabilityBounds, ~, audio2nlg] = getCellInfo(vd,cell_k,'stabilityBounds');
                shared_nlg_pulse_times = audio2nlg.shared_nlg_pulse_times - audio2nlg.first_nlg_pulse_time;
                
                stabilityBounds = 1e-3*stabilityBounds;
                sessionBounds = 1e-3*shared_nlg_pulse_times([1 end]);
                
                sessionBounds = [max(stabilityBounds(1),sessionBounds(1)) min(stabilityBounds(2),sessionBounds(2))];
                total_session_length = abs(diff(sessionBounds));
                
            elseif ~isempty(varargin)
                cut_call_data = varargin{1};
                stabilityBounds = getCellInfo(vd,cell_k,'stabilityBounds');
                
                stabilityBounds = 1e-3*stabilityBounds;
                sessionBounds = 1e-3*vertcat(cut_call_data([1 end]).corrected_callpos);
                sessionBounds = sessionBounds([1 end]);
                
                sessionBounds = [max(stabilityBounds(1),sessionBounds(1)) min(stabilityBounds(2),sessionBounds(2))];
                total_session_length = abs(diff(sessionBounds));
                
            end
            
        end
        function [fr_to_session_calls,percent_high_fr_calls,percent_session_calls,cell_ks] = get_fr_to_session_call_ratio_manual(vd,cData,b)
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
        function [fr_to_session_calls,percent_high_fr_calls,percent_session_calls,cell_ks] = get_fr_to_session_call_ratio(vd,cData,cell_ks)
            percent_session_calls = nan(1,length(cell_ks));
            percent_high_fr_calls = nan(1,length(cell_ks));
            for k = 1:length(cell_ks)
                event_pos_data = get_event_pos_high_fr(vd,cell_ks(k));
                if ~isempty(event_pos_data)
                    
                    [~, cut_call_data] = getCellInfo(vd,cell_ks(k),'cut_call_data');
                    
                    total_call_length = get_total_call_length(vd,cData,cell_ks(k),2);
                    total_session_length = get_total_session_length(vd,cell_ks(k));
                    if total_session_length < 60*60
                        total_session_length_tmp = get_total_session_length(vd,cell_ks(k),cut_call_data);
                        if total_session_length_tmp > total_session_length
                            total_session_length = total_session_length_tmp;
                        end
                    end
                    
                    percent_session_calls(k) = 100*total_call_length/total_session_length;
                    
                    cut_call_data = cut_call_data(~[cut_call_data.noise]);
                    callPos = vertcat(cut_call_data.corrected_callpos);
                    isCall = false(1,length(event_pos_data));
                    for event_k = 1:length(event_pos_data)
                        high_fr_eventpos = event_pos_data(event_k).corrected_eventpos + [-1e3 1e3];
                        intersects = (high_fr_eventpos(2) > callPos(:,1)) & (high_fr_eventpos(1) < callPos(:,2));
                        isCall(event_k) = any(intersects);
                    end
                    percent_high_fr_calls(k) = 100*sum(isCall)/length(isCall);
                end
            end
            
            fr_to_session_calls = percent_high_fr_calls./percent_session_calls;
            fr_to_session_calls(fr_to_session_calls==0) = min(fr_to_session_calls);
        end
        function event_pos_data = get_event_pos_high_fr(vd,cell_k)
            
            timestamps = getSpikes(vd,cell_k);
            timestamps = 1e-3*timestamps;
            baselineMultiplier = 0;
            min_peak_separation_s = 2;
            nPeaks = 250;
            fr_dT = 0.1;
            fs = 1/fr_dT;
            smoothSpan = (2*min_peak_separation_s)/fr_dT;
            
            [~, sessionBounds] = get_total_session_length(vd,cell_k);
            if range(timestamps)/range(sessionBounds) > 5
                [~, cut_call_data] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
                [~, sessionBounds] = get_total_session_length(vd,cell_k,cut_call_data);
            end
            timestamps = inRange(timestamps,sessionBounds);
            if isempty(timestamps)
                event_pos_data = [];
                return
            end
            t = min(timestamps):fr_dT:max(timestamps);
            
            binned_spikes = histcounts(timestamps,[t t(end)+1])/fr_dT;
            fr = smoothdata(binned_spikes,'gaussian',smoothSpan);
            
            [pks, locs] = findpeaks(fr,'MinPeakHeight',vd.avgBaseline(cell_k) + baselineMultiplier*vd.devBaseline(cell_k),'MinPeakDistance',min_peak_separation_s*fs);
            if length(pks)<nPeaks
                event_pos_data = [];
                return
            end
            [~,idx] = sort(pks);
            topIdx = ismember(locs,locs(idx(end-nPeaks:end)));
            locs = locs(topIdx);
            
            b = vd.batIdx(cell_k);
            if strcmp(vd.expType{b},'juvenile')
                ttStr = 'TT';
                audio_dir = fullfile(vd.baseDirs{b},['bat' vd.batNum{cell_k}],['neurologger_recording' vd.cellInfo{cell_k}(1:strfind(vd.cellInfo{cell_k},ttStr)-1)],'audio','ch1');
            elseif any(strcmp(vd.expType{b},{'adult','adult_operant'}))
                baseDir = fullfile('Y:\users\maimon\', [vd.expType{b} '_recording\']);
                switch vd.callType
                    case 'call'
                        audio_dir = fullfile(baseDir,datestr(vd.expDay(cell_k),'mmddyyyy'),'audio','communication','ch1');
                    case 'operant'
                        audio_dir = fullfile(baseDir,datestr(vd.expDay(cell_k),'mmddyyyy'),'operant',['box' vd.boxNums{b}]);
                end
            end
            high_activity_t = 1e3*t(locs);
            event_pos_data = struct('file_event_pos',[],'cut',[],'corrected_eventpos',[],'f_num',[],'fName',[],'fs',[],'noise',[],'expDay',[]);
            
            if strcmp(vd.callType,'operant')
                fs = 192e3;
            else
                fs = 250e3;
            end
            
            [~, ~, audio2nlg] = getCellInfo(vd,cell_k);
            nlg_offset = audio2nlg.first_nlg_pulse_time - audio2nlg.shared_nlg_pulse_times(1);
            
            for t_k = 1:length(high_activity_t)
                nlg_time = high_activity_t(t_k);
                [fName, f_num, file_event_pos] = get_avi_time_from_nlg(audio2nlg,audio_dir,nlg_time + nlg_offset,vd.callType);
                
                if ~isempty(fName)
                    
                    event_pos_data(t_k).file_event_pos = file_event_pos;
                    event_pos_data(t_k).cut = [];
                    event_pos_data(t_k).corrected_eventpos = repmat(nlg_time,1,2);
                    event_pos_data(t_k).f_num = f_num;
                    event_pos_data(t_k).fName = fName;
                    event_pos_data(t_k).fs = fs;
                    event_pos_data(t_k).noise = false;
                    event_pos_data(t_k).expDay = vd.expDay(cell_k);
                end
            end
            
        end
        
        function [r,R] = trial_fr_consistency(vdList,cell_ks,varargin)
            
            pnames = {'timeLim','nBoot','minCall'};
            dflts  = {vdList(1).spikeRange,1e3,vdList(1).minCalls};
            [timeLims,nSubsample,minCall] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            [~,t_idx] = inRange(vdList(1).time,timeLims);
            [r,R] = deal(nan(1,length(cell_ks)));
            k = 1;
            for cell_k = cell_ks
                rSubsample = nan(1,nSubsample);
                if length(vdList) == 1
                    trialFR = vdList.trialFR(cell_k);
                    used_calls = vdList.usedCalls{cell_k};
                    n_used_calls = sum(used_calls);
                    if n_used_calls < minCall
                        k = k + 1;
                        continue
                    end
                    trialFR = trialFR(used_calls,t_idx)';
                    
                    trainIdx = 1:2:minCall;
                    testIdx = setdiff(1:minCall,trainIdx);
                    
                    for boot_k = 1:nSubsample
                        callIdx = randperm(n_used_calls,minCall);
                        sample_trial_fr = trialFR(:,callIdx);
                        rSubsample(boot_k) = corr(mean(sample_trial_fr(:,trainIdx),2),mean(sample_trial_fr(:,testIdx),2));
                    end
                    
                    pairwiseR = corr(trialFR);
                    mask = tril(true(size(pairwiseR)),-1);
                    pairwiseR = pairwiseR(mask);
                   
                elseif length(vdList) == 2
                    [trialFRs,used_calls] = deal(cell(1,2));
                    for vd_k = 1:2
                        trialFR = vdList(vd_k).trialFR(cell_k{1}(vd_k));
                        used_calls{vd_k} = vdList(vd_k).usedCalls{cell_k{1}(vd_k)};
                        trialFRs{vd_k} = trialFR(used_calls{vd_k},t_idx)';
                    end
                    
                    if any(cellfun(@sum,used_calls) < minCall)
                        k = k + 1;
                        continue
                    end
                    
                    for boot_k = 1:nSubsample
                        sample_trial_fr = cell(1,2);
                        for vd_k = 1:2
                            callIdx = randperm(sum(used_calls{vd_k}),minCall);
                            sample_trial_fr{vd_k} = trialFRs{vd_k}(:,callIdx);
                        end
                        rSubsample(boot_k) = corr(mean(sample_trial_fr{1},2),mean(sample_trial_fr{2},2));
                    end
                    
                    pairwiseR = corr(trialFRs{:});
                    pairwiseR = pairwiseR(:);
                    
                end

                r(k) = nanmean(rSubsample);
                R(k) = mean(pairwiseR(~isinf(pairwiseR) & ~isnan(pairwiseR)));
                
                k = k + 1;
            end
        end
        
        function timestamps = getSpikes(vd,cell_k)
            b = vd.batIdx(cell_k);
            [stabilityBounds, ~, audio2nlg, ttDir, spikeDir] = getCellInfo(vd,cell_k,'stabilityBounds');
            
            if any(strcmp(vd.expType{b},{'adult','adult_operant'}))
                try
                    
                    timestamps = csvread(spikeDir);
                    
                catch err
                    
                    if strcmp(err.identifier,'MATLAB:textscan:EmptyFormatString')
                        timestamps = [];
                    else
                        rethrow(err)
                    end
                    
                end
                
                if size(timestamps,1) ~= 1
                    timestamps = timestamps';
                end
                
            elseif strcmp(vd.expType{b},'adult_wujie')
                
                timestamps = csvread([vd.spike_data_dir{b} vd.batNum{cell_k} '_' vd.cellInfo{cell_k} '.csv']);
                if size(timestamps,1) ~= 1
                    timestamps = timestamps';
                end
                
            elseif strcmp(vd.expType{b},'juvenile')
                timestamps = Nlx2MatSpike(ttDir,[1 0 0 0 0],0,1,[]); % load sorted cell data
                timestamps = 1e-3*timestamps - audio2nlg.first_nlg_pulse_time; % convert to ms and align to first TTL pulse on the NLG
                timestamps = inRange(timestamps,stabilityBounds);
            end
        end
        function [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = get_cell_info(vd,cell_k)
            [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
        end
        function trialFR = trialFR(vd,cell_k)
            switch vd.kernelType
                case 'manual'
                    
                    bw = vd.frBandwidth{cell_k};
                    f = gaussFilter(vd.time,bw);
                    spike_train = full(vd.trial_spike_train{cell_k});
                    trialFR = conv2(1,f,spike_train,'same')/vd.dT;
                    
                case 'fixedOpt'
                    
                    trialFR = vd.stored_trial_fr{cell_k};
                    
                case 'varOpt'
                    
                    trialFR = vd.stored_trial_fr{cell_k};
                    
            end
            
        end
        
        function batIdx = batIdx(vd,cell_k)
            batIdx = strcmp(vd.batNums,vd.batNum{cell_k});
        end
        function responsiveCells = get.responsiveCells(obj)
            responsiveCells = find(~isnan(obj.latency));
        end
        function responsiveIdx = get.responsiveIdx(obj)
            responsiveIdx = ~isnan(obj.latency);
        end
        function responsive_cells_by_bat = get.responsive_cells_by_bat(obj)
            responsive_cells_by_bat = cell(1,length(obj.batNums));
            for n = 1:length(obj.batNums)
                responsive_cells_by_bat{n} = intersect(find(strcmp(obj.batNums{n},obj.batNum)),obj.responsiveCells);
            end
        end
        function time = get.time(obj)
            time = obj.callRange(1):obj.dT:obj.callRange(2);
        end
        function nCells = get.nCells(obj)
           nCells = length(obj.cellInfo); 
        end
        function cellTable = get.cellTable(obj)
            cellTable = table(obj.cellInfo',obj.batNum','VariableNames',{'cellInfo','batNum'});
        end
        
        function vd = update_usedCalls_by_callData(vd,cData,selecFun)
            
            for cell_k = find(vd.usable)
                try
                    [stabilityBounds, cut_call_data] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
                    callNums = vd.callNum{cell_k};
                    selectedCalls = selecFun(cData,callNums)';
                    [~,~,~,~,usedCallNums] = get_used_call_spikes(vd,stabilityBounds,cut_call_data,'selectedCalls',selectedCalls);
                    vd.usedCalls{cell_k} = usedCallNums;
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
                
                switch vd.sortingMetric
                    case 'sortingQuality'
                        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
                    case 'ratio_and_distance'
                        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
                end
                if wellSorted
                    timestamps = getSpikes(vd,cell_k);
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
        function [idx,propertyNames] = compare_vd_properties(vd,vd2)
            propertyNames = properties(vd);
            assert(all(strcmp(propertyNames,properties(vd2))))
            
            metaVD = metaclass(vd);
            dependentIdx = arrayfun(@(x) x.HasDefault,metaVD.PropertyList);
            
            propertyNames = propertyNames(dependentIdx);
            nProperty = length(propertyNames);
            idx = false(1,nProperty);
            for p_k = 1:length(propertyNames)
                switch class(vd.(propertyNames{p_k}))
                    
                    case 'double'
                        
                        idx(p_k) = all(vd.(propertyNames{p_k}) == vd2.(propertyNames{p_k}));
                        
                    case 'char'
                        
                        idx(p_k) = all(strcmp(vd.(propertyNames{p_k}),vd2.(propertyNames{p_k})));
                        
                    case 'logical'
                        
                        idx(p_k) = all(vd.(propertyNames{p_k}) == vd2.(propertyNames{p_k}));
                        
                end
            end
            
        end
        function vd = join(vd1,vd2)
            idx = compare_vd_properties(vd1,vd2);
            assert(all(idx));
            
            vd = vd1;
            propertyNames = properties(vd1);
            
            metaVD = metaclass(vd);
            dependentIdx = arrayfun(@(x) x.HasDefault || x.Dependent || strcmp(x.SetAccess,'private'),metaVD.PropertyList);
            propertyNames = propertyNames(~dependentIdx);
            
            superclass_property_names = arrayfun(@(prop) prop.Name,metaVD.SuperclassList.PropertyList,'un',0);
            propertyNames = propertyNames(~ismember(propertyNames,superclass_property_names));
            
            propertyLength = cellfun(@(x) length(vd1.(x)),propertyNames);
            propertyNames  = propertyNames(propertyLength == vd1.nCells);
            
            n_prop_to_join = length(propertyNames);
            for prop_k = 1:n_prop_to_join
                pName = propertyNames{prop_k};
                sz = size(vd1.(pName));
                joinDim = find(sz == vd1.nCells);
                vd.(pName) = cat(joinDim,vd.(pName),vd2.(pName));
            end           
            
            nBat = length(vd1.batNums);
            superclass_property_length = cellfun(@(x) length(vd1.(x)),superclass_property_names);
            str_property_idx = arrayfun(@(prop) ischar(vd1.(prop.Name)),metaVD.SuperclassList.PropertyList);
            bat_superclass_property_idx = superclass_property_length == nBat & ~str_property_idx;
            
            propIdx = find(bat_superclass_property_idx);
            for prop_k = 1:length(propIdx)
                pName = superclass_property_names{propIdx(prop_k)};
                sz = size(vd1.(pName));
                joinDim = find(sz == nBat);
                vd.(pName) = cat(joinDim,vd.(pName),vd2.(pName));
            end
            
            propIdx = find(~bat_superclass_property_idx);
            for prop_k = 1:length(propIdx)
                pName = superclass_property_names{propIdx(prop_k)};
                if ~all(size(vd1.(pName)) == size(vd2.(pName))) || ~isequal(vd1.(pName),vd2.(pName))
                   vd.(pName) = []; 
                end
            end
           

        end
        function [cellTable,idx1,idx2] = cellIntersect(vd1,vd2)
            vdList = [vd1 vd2];
            usableIdx = cell(1,2);
            for vd_k = 1:2
                usableIdx{vd_k} = find(vdList(vd_k).usable);
            end
            [cellTable,idx1,idx2] = intersect(vd1.cellTable(vd1.usable,:),vd2.cellTable(vd2.usable,:));
            idx1 = usableIdx{1}(idx1);
            idx2 = usableIdx{2}(idx2);
        end
        function [shared_resp_cell_table,idx1,idx2] = shared_responsive_cells(vd1,vd2)
            shared_resp_cell_table = cell(1,2);
            [shared_cell_table,idx1,idx2] = cellIntersect(vd1,vd2);
            vdList = [vd1 vd2];
            idxList = {idx1, idx2};
            for vd_k = 1:2
                resp_idx = vdList(vd_k).responsiveIdx(idxList{vd_k});
                shared_resp_cell_table{vd_k} = shared_cell_table(resp_idx,:);
                idxList{vd_k} = idxList{vd_k}(resp_idx);
            end
            [idx1,idx2] = deal(idxList{:});
        end
        
        function export_all_call_spike_data(vd,output_dir)
            lastProgress = 0;
            usedCells = false(1,vd.nCells);
            for cell_k = 1:vd.nCells
                
                switch vd.sortingMetric
                    case 'sortingQuality'
                        wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
                    case 'ratio_and_distance'
                        wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
                end
                if wellSorted
                    timestamps = getSpikes(vd,cell_k);
                    [stabilityBounds, cut_call_data] = getCellInfo(vd,cell_k,'stabilityBounds','cut_call_data');
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

function [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = getCellInfo(vd,cell_k,varargin)

b = vd.batIdx(cell_k);
cellInfo = vd.cellInfo{cell_k};
baseDir = vd.baseDirs{b};
batNum = vd.batNums{b};

if any(strcmp(vd.expType{b},{'adult','adult_operant'}))
    
    spikeDir = fullfile(baseDir,'spike_data',[batNum '_' cellInfo '.csv']);
    %       If we want to calculate spike waveform stats, need to include path to actual .ntt file
    ttDir = [];
    call_data_dir = fullfile(baseDir,'call_data');
    exp_date_str = cellInfo(1:8);
    
    try
        if any(strcmp(vd.callType,{'call','echo'}))
            audio2nlg = load(fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit.mat'])); % load fit data to sync audio to nlg data
        elseif strcmp(vd.callType,'operant')
            s = load(fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit.mat']),'first_nlg_pulse_time','first_audio_pulse_time');
            audio2nlg = load(fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit_operant_box_'  vd.boxNums{b} '.mat'])); % load fit data to sync audio to nlg data
            audio2nlg.first_nlg_pulse_time = s.first_nlg_pulse_time;
            audio2nlg.first_audio_pulse_time = s.first_audio_pulse_time;
        end
    catch
        [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = deal([]);
        return
    end
    
    if any(strcmp(varargin,'stabilityBounds')) && strcmp(vd.cellType,'singleUnit')
        s = load(fullfile(vd.analysisDir{b},'cell_stability_info.mat'));
        cell_stability_info = s.cell_stability_info;
        idx = find(strcmp(cellInfo,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
        stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
    else
        stabilityBounds = [];
    end
    
    if any(strcmp(varargin,'cut_call_data'))
        
        if any(strcmp(vd.callType,{'call','operant'}))
            
            if strcmp(vd.callType,'call')
                cut_call_fname = fullfile(call_data_dir,[exp_date_str '_cut_call_data.mat']);
            elseif strcmp(vd.callType,'operant')
                cut_call_fname = fullfile(call_data_dir,[exp_date_str '_cut_call_data_operant_box_' vd.boxNums{b} '.mat']);
                if any(strcmp(vd.operant_reward_status,{'rewardedOnly','notRewarded'}))
                    operant_reward_delay = 5e3;
                    operantEvents = get_operant_events(vd,cell_k);
                end
            end
            
            if exist(cut_call_fname,'file')
                s = load(cut_call_fname);
                cut_call_data = s.cut_call_data;
            else
                [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = deal([]);
                return
            end
            
            if strcmp(vd.selectCalls,'selfCall')
                used_call_idx =  find(strcmp({cut_call_data.batNum},batNum));
                not_used_call_idx = find(~strcmp({cut_call_data.batNum},batNum) & ~strcmp({cut_call_data.batNum},'noise'));
            elseif  strcmp(vd.selectCalls,'otherCall')
                used_call_idx = find(~strcmp({cut_call_data.batNum},batNum) & ~strcmp({cut_call_data.batNum},'noise'));
                not_used_call_idx = find(strcmp({cut_call_data.batNum},batNum));
            elseif strcmp(vd.selectCalls,'allCall')
                used_call_idx = find(~strcmp({cut_call_data.batNum},'noise'));
                not_used_call_idx = [];
            end
            
            if vd.exclude_neighboring_calls && ~isempty(used_call_idx) && ~isempty(not_used_call_idx)
                used_call_pos = vertcat(cut_call_data(used_call_idx).corrected_callpos);
                not_used_bat_call_pos = vertcat(cut_call_data(not_used_call_idx).corrected_callpos) + 1e3*vd.spikeRange;
                not_used_bat_call_pos = reshape(not_used_bat_call_pos,[],1);
                
                nCall = length(used_call_idx);
                neighboring_call_idx  = false(1,nCall);
                
                for call_k = 1:nCall
                    cp = used_call_pos(call_k,:);
                    neighboring_call_idx(call_k) = any(diff(sign(repmat(cp,size(not_used_bat_call_pos,1),1) - repmat(not_used_bat_call_pos,1,2)),[],2));
                end
                used_call_idx = used_call_idx(~neighboring_call_idx);
            end
            
            if strcmp(vd.callType,'operant') && any(strcmp(vd.operant_reward_status,{'rewardedOnly','notRewarded'}))
                if ~isempty(operantEvents)
                    rewarded_event_idx = operantEvents.food_port_status.FoodPortBack | operantEvents.food_port_status.FoodPortFront;
                    rewarded_event_times = operantEvents.eventTimes(rewarded_event_idx);
                    used_call_pos = vertcat(cut_call_data(used_call_idx).corrected_callpos);
                    
                    nCall = length(used_call_idx);
                    dT_call = nan(1,nCall);
                    
                    for call_k = 1:nCall
                        dT_call(call_k) = min(abs(used_call_pos(call_k,1) - rewarded_event_times));
                    end
                    if strcmp(vd.operant_reward_status,'rewardedOnly')
                        used_call_idx = used_call_idx(dT_call < operant_reward_delay);
                    elseif strcmp(vd.operant_reward_status,'notRewarded')
                        used_call_idx = used_call_idx(dT_call > operant_reward_delay);
                    end
                else
                    used_call_idx = [];
                end
            end
            
            cut_call_data = cut_call_data(used_call_idx);
            
        elseif strcmp(vd.callType,'echo')
            cut_echo_fname = fullfile(call_data_dir,[exp_date_str '_cut_echo_data.mat']);
            if exist(cut_echo_fname,'file')
                s = load(cut_echo_fname);
                cut_call_data = s.cut_call_data;
            else
                [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = deal([]);
                return
            end
        elseif strcmp(vd.callType,'operant_reward')
            
            operantEvents = get_operant_events(vd,cell_k);
            
            if ~isempty(operantEvents)
                
                rewarded_event_idx = (operantEvents.food_port_status.FoodPortBack | operantEvents.food_port_status.FoodPortFront) & ~isnan(operantEvents.delay2reward) &  ~isnan(operantEvents.eventTimes);
                nReward = sum(rewarded_event_idx);
                rewardTimes = operantEvents.eventTimes(rewarded_event_idx) + 1e3*operantEvents.delay2reward(rewarded_event_idx);
                rewardTimes = repmat(rewardTimes,1,2);
                cut_call_data = struct('corrected_callpos',mat2cell(rewardTimes,ones(1,nReward)),'uniqueID',num2cell(-nReward:-1)','batNum',repmat({batNum},nReward,1));
                
            else
                cut_call_data = [];
            end
            
            
        end
    else
        cut_call_data = [];
    end
    
    stabilityBounds = stabilityBounds - audio2nlg.first_nlg_pulse_time;
    success = true;
    
elseif strcmp(vd.expType{b},'adult_wujie')
    
    expDate = cellInfo(1:strfind(cellInfo,vd.tetrodeStr)-1);
    audioDir = [baseDir 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
    stabilityBounds = [-Inf Inf];
    
    try
        if any(strcmp(varargin,'cut_call_data'))
            s = load([audioDir 'cut_call_data.mat']);
            cut_call_data = s.cut_call_data;
            cut_call_data = cut_call_data(~[cut_call_data.noise]);
            if isempty(cut_call_data)
                audio2nlg = [];
                ttDir = [];
                success = true;
                return
            end
            callIdx = unique(cellfun(@(call) find(cellfun(@(bNum) strcmp(bNum,vd.batNums{b}),call)),{cut_call_data.batNum}));
            
            if length(callIdx) == 1
                callpos = horzcat(cut_call_data.corrected_callpos);
                callpos = callpos(callIdx,:);
                [cut_call_data.corrected_callpos] = deal(callpos{:});
            else
                keyboard
            end
            call_info_fname = dir([audioDir 'call_info_*_' vd.callType '_' expDate '.mat']);
            if length(call_info_fname) > 1
                keyboard
            end
            s = load(fullfile(call_info_fname.folder,call_info_fname.name));
            call_info = s.call_info;
            
            assert(all([cut_call_data.uniqueID] == [call_info.callID]));
            
            bat_calls = cellfun(@(x) ischar(x{1}) && contains(x,batNum),{call_info.behaviors});
            cut_call_data = cut_call_data(bat_calls);
        else
            cut_call_data = [];
        end
        
        audio2nlg = [];
        ttDir = [];
        spikeDir = [];
        
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
    
elseif strcmp(vd.expType{b},'juvenile')
    
    cluster_num_idx = strfind(cellInfo,vd.clusterStr)+length(vd.clusterStr);
    tt_num_idx = strfind(cellInfo,vd.tetrodeStr )+length(vd.tetrodeStr );
    
    expDate = cellInfo(1:strfind(cellInfo,vd.tetrodeStr)-1);
    tetrodeNum = str2double(cellInfo(tt_num_idx));
    cellNum = str2double(cellInfo(cluster_num_idx:end));
    
    if cellNum < 10
        sorted_cell_string = [vd.clusterStr '0'];
    else
        sorted_cell_string = vd.clusterStr;
    end
    
    tetrode = [expDate vd.tetrodeStr num2str(tetrodeNum) sorted_cell_string num2str(cellNum)]; % build tetrode/cell string
    ttDir = [vd.spike_data_dir{b} 'bat' batNum filesep expDate filesep tetrode '.ntt']; % filename for cell of interest
    audioDir = [baseDir 'bat' batNum filesep 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
    
    audio2nlg = load([audioDir 'audio2nlg_fit.mat']); % load fit data to sync audio to nlg data
    
    if any(strcmp(varargin,'stabilityBounds'))
        s = load(fullfile(vd.analysisDir{b},'cell_stability_info.mat'));
        cell_stability_info = s.cell_stability_info;
        idx = find(strcmp(cellInfo,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
        stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
        stabilityBounds = stabilityBounds - audio2nlg.first_nlg_pulse_time;
    else
        stabilityBounds = [];
    end
    
    
    try
        if any(strcmp(varargin,'cut_call_data'))
            if strcmp(vd.callType,'call') % are we looking at calls or echolocation clicks?
                s = load([audioDir 'cut_call_data.mat']);
                cut_call_data = s.cut_call_data;
                cut_call_data = cut_call_data(~[cut_call_data.noise]);
                
            elseif strcmp(vd.callType,'echo')
                s = load([audioDir 'cut_echo_data.mat']);
                cut_call_data = s.cut_call_data;
                cut_call_data = cut_call_data(~[cut_call_data.noise]);
                
                s = load([audioDir 'juv_call_info_' vd.batNums{b} '_echo.mat']);
                call_info = s.call_info;
                
                assert(all([cut_call_data.uniqueID] == [call_info.callID]));
                
                echo_calls = strcmp({call_info.echoCall},'juvEcho');
                cut_call_data = cut_call_data(echo_calls);
            end
        else
            cut_call_data = [];
        end
        
    catch err
        
        disp(err)
        keyboard
        
        success = input('continue?');
        
        if ~success
            return
        end
        
    end
    
    success = true;
    spikeDir = [];
end

end

function operantEvents = get_operant_events(vd,cell_k)

b = vd.batIdx(cell_k);
cellInfo = vd.cellInfo{cell_k};
baseDir = vd.baseDirs{b};

call_data_dir = fullfile(baseDir,'call_data');
exp_date_str = cellInfo(1:8);

event_time_fname = fullfile(call_data_dir,[exp_date_str '_event_times_box_' vd.boxNums{b} '.mat']);
if exist(event_time_fname,'file')
    operantEvents = load(event_time_fname);
else
    operantEvents = [];
    return
end

[operantEvents.eventTimes,idx] = sort(operantEvents.eventTimes);
operantEvents.food_port_status = operantEvents.food_port_status(idx,:);
operantEvents.delay2reward = operantEvents.delay2reward(idx,:);

end

function [callSpikes, callNum, call_bat_num, callLength, usedCalls] = get_used_call_spikes(vd,stabilityBounds,cut_call_data,varargin)

pnames = {'selectedCalls','timestamps'};
dflts  = {[cut_call_data.uniqueID],[]};
[selectedCalls,timestamps] = internal.stats.parseArgs(pnames,dflts,varargin{:});

all_call_idx = ismember([cut_call_data.uniqueID],selectedCalls);
all_call_times = [cut_call_data(all_call_idx).corrected_callpos];

spikeRange = 1e3*vd.spikeRange;

nCalls = length(cut_call_data); % total # of calls, excluding any errors by cutting procedure

call_k = 1; % counter for all used calls

callSpikes = cell(1,nCalls); % initialize cell to contain spike times relative to call onset
callNum = nan(1,nCalls);
call_bat_num = cell(1,nCalls);
callLength = nan(1,nCalls);
usedCalls = true(1,nCalls);

for call = 1:nCalls % iterate through all the calls within this .wav file
    cp = cut_call_data(call).corrected_callpos; % time of onset and offset of the call of interest on this iteration
    callposOffset = [cp(1)+spikeRange(1), cp(2)+spikeRange(2)]; % total time window around call of interest
    
    if ~any(isnan(cp)) % check that the TTL alignment worked for this call
        
        in_stable_range = ~strcmp(vd.cellType,'singleUnit'); % Initialize as false only if we're using single units (MUA is always "stable")
        for bound_k = 1:size(stabilityBounds,1) % check that this call is within timeframe that this cell is stable
            in_stable_range = in_stable_range || (callposOffset(1) > stabilityBounds(bound_k,1) && callposOffset(2) < stabilityBounds(bound_k,2));
        end
        
        if in_stable_range
            
            if (strcmp(vd.onset_or_offset,'onset') && any(all_call_times > (cp(1) + spikeRange(1)) & all_call_times < cp(1))) || ... skip any calls within a call train
                    (strcmp(vd.onset_or_offset,'offset') && any(all_call_times < (cp(2) + spikeRange(2)) & all_call_times > cp(2))) || ...
                    ~(ismember(cut_call_data(call).uniqueID,selectedCalls))
                
                usedCalls(call_k) = false;
            end
            
            if ~isempty(timestamps)
                relativeCP = cp(strcmp({'onset','offset'},vd.onset_or_offset)); % choose beginning or end of call
                timestampsCall = inRange(timestamps,callposOffset); % all spikes occuring within that time window
                callSpikes{call_k} = 1e-3*(timestampsCall - relativeCP); % align spikes to call onset/offset, convert to sec, and store
                callNum(call_k) = cut_call_data(call).uniqueID;
                call_bat_num{cell_k} = cut_call_data(call).batNum;
                callLength(call_k) = 1e-3*diff(cp);
                if vd.timeWarp
                    callSpikes{call_k} = callSpikes{call_k} * vd.warp_call_length/callLength(call_k);
                end
            end
            call_k = call_k + 1; % increment total calls
        end
        
    end
end

%now truncate the list of calls to how many calls actually occurred
nCalls = call_k - 1; % determine the actual number of calls
callSpikes = callSpikes(1:nCalls); % truncate list of spikes relative to call onset to length of n_calls
callNum = callNum(1:nCalls);
call_bat_num = call_bat_num(1:nCalls);
callLength = callLength(1:nCalls);
usedCalls = usedCalls(1:nCalls);
end

function [baselineMean, baselineStd] = calculate_baseline(vd,cell_k,timestamps,cut_call_data)

b = vd.batIdx(cell_k);
switch vd.baselineMethod
    case 'preCall'
        baselineLength = abs(diff(vd.baselineRange));
        %         nTrial = sum(vd.usedCalls{cell_k});
        %         allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        %         usedSpikes = inRange(allSpikes,vd.baselineRange);
        %         baselineMean = length(usedSpikes)/nTrial/baselineLength;
        baselineMean = mean(cellfun(@(x) length(inRange(x,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k})));
        baselineStd = std(cellfun(@(x) length(inRange(x,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(vd.usedCalls{cell_k})));
    case 'randomSamples'
        if nargin < 3
            [~, cut_call_data, ~, ~] = getCellInfo(vd,cell_k,'cut_call_data');
            timestamps = getSpikes(vd,cell_k);
        end
        
        if any(strcmp(vd.expType{b},{'adult','adult_operant'}))
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
            [~, cut_call_data] = getCellInfo(vd,cell_k,'cut_call_data');
            timestamps = getSpikes(vd,cell_k);
        end
        
        callRange = [cut_call_data(1).corrected_callpos(1) cut_call_data(end).corrected_callpos(2)];
        timestamps = inRange(timestamps,callRange);
        
        if any(strcmp(vd.expType{b},{'adult','adult_operant'}))
            timestamps = timestamps - timestamps(1);
        end
        
        timestamps = timestamps(timestamps>0);
        
        t = linspace(0,max(timestamps),1e4);
        dT = mean(diff(t));
        binned_spikes = histcounts(timestamps,[t t(end)+1])/(1e-3*dT);
        fr = smoothdata(binned_spikes,2,'gaussian',1e2);
        
        call_t = true(1,length(t));
        baselineLength = vd.baseline_sample_length*1e3;
        
        for k = 1:length(cut_call_data)
            call_t(t > (cut_call_data(k).corrected_callpos(1) - baselineLength) & t < (cut_call_data(k).corrected_callpos(2) + baselineLength)) = false;
        end
        
        baselineMean = mean(fr(call_t));
        baselineStd = std(fr(call_t));
end

end

function usable = checkUsability(vd,cell_k)

allSpikes = vd.callSpikes{cell_k}(vd.usedCalls{cell_k});
if isempty(allSpikes)
    usable = false;
    return
end
allSpikes = [allSpikes{:}];
allSpikes = inRange(allSpikes,vd.baselineRange);
nTrial = sum(vd.usedCalls{cell_k});

switch vd.cellType
    
    case 'singleUnit'
        
        switch vd.sortingMetric
            case 'sortingQuality'
                wellSorted = vd.sortingQuality(cell_k) <= vd.sortingThreshold;
            case 'ratio_and_distance'
                wellSorted = vd.LRatio(cell_k) <= vd.sortingThreshold(1) & vd.isolationDistance(cell_k) >= vd.sortingThreshold(2);
        end
        usable = (length(allSpikes)/sum(vd.usedCalls{cell_k}))/range(vd.baselineRange) >= vd.minSpikes &&...
            nTrial >= vd.minCalls &&...
            wellSorted;
        
    case 'multiUnit'
        
         usable = (length(allSpikes)/sum(vd.usedCalls{cell_k}))/range(vd.baselineRange) >= vd.minSpikes && nTrial >= vd.minCalls;
        
end

end

function [spikeRate, spikeTrain, bw] = frKernelEstimate(vd,cell_k)

nTrial = length(vd.callSpikes{cell_k});
nT = length(vd.time);
spikeTrain = zeros(nTrial,nT);


switch vd.kernelType
    case 'manual'
        
        if isnan(vd.constantBW)
            bw = kde([vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}]);
        else
            bw = vd.constantBW;
        end
        
        for tt = 1:nTrial
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[vd.time(1),vd.time(end)]);
            if ~isempty(usedSpikes)
                spikeTrain(tt,:) = histcounts(usedSpikes,[vd.time vd.time(end)+vd.dT]);
            end
        end
        
        spikeRate = [];
        
    case 'fixedOpt'
        spikeRate = zeros(nTrial,nT);
        bw = nan(nTrial,1);
        for tt = 1:nTrial
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[vd.time(1),vd.time(end)]);
            if ~isempty(usedSpikes)
                spikeTrain(tt,:) = histcounts(usedSpikes,[vd.time vd.time(end)+vd.dT]);
                [spikeRate(tt,:), ~, bw(tt)]= sskernel(usedSpikes,vd.time,vd.frBandwidthRange);
            end
        end
        
    case 'varOpt'
        spikeRate = zeros(nTrial,nT);
        bw = nan(nTrial,nT);
        for tt = 1:nTrial
            usedSpikes = [vd.callSpikes{cell_k}{tt}];
            usedSpikes = inRange(usedSpikes,[vd.time(1),vd.time(end)]);
            if ~isempty(usedSpikes)
                spikeTrain(tt,:) = histcounts(usedSpikes,[vd.time vd.time(end)+vd.dT]);
                [spikeRate(tt,:), ~, bw(tt,:)]= ssvkernel(usedSpikes,time,vd.frBandwidthRange);
            end
        end
        
        
end

spikeTrain = sparse(spikeTrain);

end

function [spike_rate_mean, spike_rate_dev] = calculateAvgFR(vd,cell_k)

switch vd.kernelType
    case 'manual'
        
        spikeRate = vd.trialFR(cell_k);
        spikeRate = spikeRate(vd.usedCalls{cell_k},:);
        
        spike_rate_mean = mean(spikeRate);
        spike_rate_dev = std(spikeRate)/sqrt(size(spikeRate,1))';
        spike_rate_dev = repmat(spike_rate_dev,2,1);
        
    case 'fixedOpt'
        
        allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        [spike_rate_mean,~,~,~,~,spike_rate_dev] = sskernel(allSpikes,vd.time);
        spike_rate_mean = spike_rate_mean/sum(vd.usedCalls{cell_k});
        
    case 'varOpt'
        
        allSpikes = [vd.callSpikes{cell_k}{vd.usedCalls{cell_k}}];
        [spike_rate_mean,~,~,~,~,spike_rate_dev] = ssvkernel(allSpikes,vd.time);
        spike_rate_mean = spike_rate_mean/sum(vd.usedCalls{cell_k});
        
end

end

function f = gaussFilter(t,sigma)

f = exp( -t.^2 / (2*sigma^2));
f = f / sum(f);

end

function [latency, respType, respValency, respStrength] = calculateLatency(vd,cell_k,varargin)

pnames = {'trialIdx','latencyType'};
dflts  = {[],vd.latencyType};
[trialIdx,latencyType] = internal.stats.parseArgs(pnames,dflts,varargin{:});


if any(strcmp(latencyType,{'pre','post','during','prePost','peri'}))
    [responsive, respType, respValency, respStrength] = periCallResponsive(vd,cell_k,trialIdx);
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
    
elseif strcmp(latencyType,'cusum')
    [respType, respValency, respStrength, latency] = cusumResponsive(vd,cell_k,trialIdx);
    if isempty(inRange(latency,vd.latencyRange))
        latency = NaN;
        respType = NaN;
        respValency = NaN;
        respStrength = NaN;
    end
    
elseif strcmp(latencyType,'trialBased')
    [respType, respValency, respStrength, latency]  = trialBasedResponsive(vd,cell_k,trialIdx);
elseif strcmp(latencyType,'std_above_baseline')
    [respType, respValency, respStrength, latency]  = std_above_baseline_responsive(vd,cell_k,trialIdx);
elseif strcmp(latencyType,'sliding_win_test')
    [respType, respValency, respStrength, latency]  = sliding_win_test_responsive(vd,cell_k,trialIdx);
end


end

function [respType, respValency, respStrength, latency] = sliding_win_test_responsive(vd,cell_k,trialIdx)

if nargin < 3 || isempty(trialIdx)
    trialIdx = vd.usedCalls{cell_k};
end

[respType, respValency, respStrength, latency] = deal(NaN);
[latency_time,time_idx] = inRange(vd.time,vd.latencyRange);
nT = sum(time_idx);
win_size_s = 0.25;
winSize = win_size_s/vd.dT;
slidingWinIdx = slidingWin(nT,winSize,winSize-1);

baselineLength = diff(vd.baselineRange);
baseSpikes = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/baselineLength,vd.callSpikes{cell_k}(trialIdx));
baseFR = mean(baseSpikes);

trialSpikes = full(vd.trial_spike_train{cell_k}(trialIdx,time_idx));
avgFR = sum(trialSpikes)/size(trialSpikes,1)/vd.dT;
slidingFR = avgFR(slidingWinIdx);

[~,p] = ttest(slidingFR', baseFR);

h = p < vd.responsiveAlpha/size(slidingWinIdx,1);

if any(h)
    [latency,respType] = get_thresh_crossing(vd,h,latency_time);
    if ~isnan(latency)
        respStrength = median(avgFR(h)) - vd.avgBaseline(cell_k);
        respValency = sign(respStrength);
    end
end



end

function [respType, respValency, respStrength, latency] = cusumResponsive(vd,cell_k,trialIdx)

if nargin < 3 || isempty(trialIdx)
    trialIdx = vd.usedCalls{cell_k};
end

spikeTrain = sum(vd.trial_spike_train{cell_k}(trialIdx,:),1);
idx = find(vd.time>=vd.cusumStart & vd.time<vd.latencyRange(2));
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

function [responsive, respType, respValency, respStrength] = periCallResponsive(vd,cell_k,trialIdx)

if nargin < 3 || isempty(trialIdx)
    trialIdx = vd.usedCalls{cell_k};
end

usedSpikes = vd.callSpikes{cell_k}(trialIdx);

baselineLength = diff(vd.baselineRange);
baseSpikes = cellfun(@(spikes) length(inRange(spikes,vd.baselineRange))/baselineLength,usedSpikes);

switch vd.latencyType
    
    case 'during'
        callSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,usedSpikes,num2cell(vd.callLength{cell_k}(usedSpikes)));
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'during';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'pre'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,usedSpikes);
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'pre';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'post'
        callSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,usedSpikes);
        responsive = ttest(baseSpikes,callSpikes,'Alpha',vd.responsiveAlpha);
        if responsive
            respValency = sign(mean(callSpikes) - mean(baseSpikes));
            respType = 'post';
            respStrength = abs(mean(callSpikes) - mean(baseSpikes));
        end
    case 'prePost'
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,usedSpikes);
        postCallSpikes = cellfun(@(spikes) length(inRange(spikes,[0,vd.postCall]))/vd.postCall,usedSpikes);
        responsivePre = ttest(baseSpikes,preCallSpikes,'Alpha',vd.responsiveAlpha);
        responsivePost = ttest(baseSpikes,postCallSpikes,'Alpha',vd.responsiveAlpha);
        responsive = (~isnan(responsivePre) && ~isnan(responsivePost)) && (responsivePre || responsivePost);
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
        preCallSpikes = cellfun(@(spikes) length(inRange(spikes,[-vd.preCall,0]))/vd.preCall,usedSpikes);
        postCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[callLength,callLength + vd.postCall]))/vd.postCall,usedSpikes,num2cell(vd.callLength{cell_k}(trialIdx)));
        duringCallSpikes = cellfun(@(spikes,callLength) length(inRange(spikes,[0,callLength]))/callLength,usedSpikes,num2cell(vd.callLength{cell_k}(trialIdx)));
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

function [respType, respValency, respStrength, latency] = trialBasedResponsive(vd,cell_k,trialIdx)

if nargin < 3 || isempty(trialIdx)
    trialIdx = vd.usedCalls{cell_k};
end

[respType, respValency, respStrength, latency] = deal(NaN);

resp_trials = get_responsive_trials(vd,cell_k,vd.responsive_nStd_over_baseline,vd.latencyRange,trialIdx);

if sum(resp_trials) >= vd.min_responsive_trials && sum(resp_trials)/length(resp_trials) > vd.min_fraction_responsive_trials
    
    [latency_time,time_idx] = inRange(vd.time,vd.latencyRange);
    
    used_calls = find(trialIdx);
    
    trialFR = vd.trialFR(cell_k);
    
    fr = mean(trialFR(used_calls(resp_trials),:));
    frDev = std(trialFR(used_calls(resp_trials),:))/sqrt(sum(resp_trials));
    
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

end

function [respType, respValency, respStrength, latency] = std_above_baseline_responsive(vd,cell_k,trialIdx)

if nargin < 3 || isempty(trialIdx)
    trialIdx = vd.usedCalls{cell_k};
end

[respType, respValency, respStrength, latency] = deal(NaN);
[latency_time,time_idx] = inRange(vd.time,vd.latencyRange);
trialFR = vd.trialFR(cell_k);

fr = mean(trialFR(trialIdx,:));
frDev = std(trialFR(trialIdx,:))/sqrt(sum(trialIdx));

thresh = vd.avgBaseline(cell_k) + vd.responsive_nStd_over_baseline*vd.devBaseline(cell_k);

latency_time_fr = fr(time_idx);

over_thresh = latency_time_fr - frDev(time_idx) > thresh;

if any(over_thresh)
    
    [latency,respType] = get_thresh_crossing(vd,over_thresh,latency_time);
    if ~isnan(latency)
        respValency = 1;
        respStrength = median(latency_time_fr(over_thresh)) - vd.avgBaseline(cell_k);
    end
end



end

function [latency,respType] = get_thresh_crossing(vd,over_thresh,latency_time)

thresh_crossing = find(diff(over_thresh) == 1) + 1;
crossing_length = zeros(1,length(thresh_crossing));

k = 1;
for tc = thresh_crossing
    if ~isempty(find(diff(over_thresh(tc:end)) == -1,1))
        crossing_length(k) = find(diff(over_thresh(tc:end)) == -1,1);
    else
        crossing_length(k) = length(vd.time) - tc;
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
else
    [latency,respType] = deal(NaN);
end

end

function [meanFR, medianISI, peak2trough, spikeWidth, avg_spike, duration] = get_spike_stats(timestamps,ttDir)
ADC_period = 34.13; % in microseconds

[~, samples] = Nlx2MatSpike(ttDir,[1 0 0 0 1],0,1,[]);

n_samples = size(samples,1);

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

height_width_intersect = repmat(halfSpikePeak,size(interp_time)) - avg_spike_interp;
spikeWidth = abs(diff(interp_time(diff(sign(height_width_intersect)) ~= 0)));

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

function merged = merge_wins(wins,fs,thresh)
if size(wins,1)>1
    diffs = (wins(2:length(wins),1)- wins(1:length(wins)-1,2));
    if any(diffs/fs < thresh)
        w = find(diffs/fs < thresh,1);
        wins(w,:) = [wins(w,1),wins(w+1,2)];
        wins_del = 1:length(wins) ~= w+1;
        merged = wins(wins_del,:);
        merged = merge_wins(merged,fs,thresh);
    else
        merged = wins;
    end
else
    merged = wins;
end

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
end

function vd = updateUsability(vd)

updated_cells = 0;
for cell_k = 1:vd.nCells
    usability = checkUsability(vd,cell_k);
    if vd.usable(cell_k) ~= usability
        if usability - vd.usable(cell_k) == 1 % cell is now usable
            [vd.stored_trial_fr{cell_k},vd.trial_spike_train{cell_k},...
                vd.frBandwidth{cell_k}] = frKernelEstimate(vd,cell_k);
            [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k);
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
disp([num2str(updated_cells) ' cells updated']);
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
    [vd.avgFR{cell_k}, vd.devFR{cell_k}] = calculateAvgFR(vd,cell_k);
    
end
vd = updateBaseline(vd);


end
