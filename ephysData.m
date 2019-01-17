classdef ephysData
    properties 
        baseDirs
        batNums 
        dateFormat
        birthDates 
        analysisDir 
        activeChannels
        spike_data_dir
        lfp_data_dir
        expType
    end
    properties
        call_echo = 'call'
    end
    
    methods
       
        function ed = ephysData(expType)
            
            if nargin == 0
                ed.expType = 'juvenile';
            else
                ed.expType = expType;
            end
            
            switch ed.expType
                
                case 'juvenile'
                    
                    ed.batNums = {'71319','71284','71173','13681'};
                    ed.baseDirs = repmat({'E:\ephys\juvenile_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = {datetime(2016,4,23),datetime(2016,09,24),datetime(2016,09,21),datetime(2017,10,20)};num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = 'C:\Users\phyllo\Documents\Maimon\ephys\data_analysis_results\';
                    ed.activeChannels = {setdiff(0:15,4),setdiff(0:15,[1:4, 8, 9]),setdiff(0:15,12),setdiff(0:15,8)};
                    ed.spike_data_dir = 'E:\ephys\juvenile_recording\tetrode_data\';
                    ed.lfp_data_dir = [];
                    
                case 'adult'
                    
                    ed.batNums = {'14620','71334','65694','71360'};
                    ed.baseDirs = repmat({'E:\ephys\adult_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = 'E:\ephys\adult_recording\data_analysis_results\';
                    ed.activeChannels = {0:15,0:15,0:15,setdiff(0:15,[7 8])};
                    ed.spike_data_dir = 'E:\ephys\adult_recording\spike_data\';
                    ed.lfp_data_dir = [];
                    
                    
                case 'adult_wujie'
                    
                    ed.batNums = {'59813','71348','65997'};
                    ed.baseDirs = repmat({'E:\ephys\adult_recording\wujie_data\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = 'E:\ephys\adult_recording\data_analysis_results\';
                    ed.activeChannels = [];
                    ed.spike_data_dir = 'E:\ephys\adult_recording\spike_data\';
                    ed.lfp_data_dir = 'E:\ephys\adult_recording\lfp_data\';
                    
            end
            
        end
        
        function timestamps = getSpikes(ed,b,cell_k)
            
            [stabilityBounds, ~, audio2nlg, ttDir, spikeDir] = getCellInfo(ed,b,cell_k);
            
            switch ed.expType
                
                case 'adult'
                    
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
                    
                case 'adult_wujie'
                    
                    timestamps = csvread([ed.spike_data_dir ed.batNum{cell_k} '_' ed.cellInfo{cell_k} '.csv']);
                    if size(timestamps,1) ~= 1
                        timestamps = timestamps';
                    end
                    
                case 'juvenile'
                    timestamps = Nlx2MatSpike(ttDir,[1 0 0 0 0],0,1,[]); % load sorted cell data
                    timestamps = 1e-3*timestamps - audio2nlg.first_nlg_pulse_time; % convert to ms and align to first TTL pulse on the NLG
                    timestamps = inRange(timestamps,stabilityBounds);
            end
        end
        
        function [stabilityBounds, cut_call_data, audio2nlg, ttDir, spikeDir, success] = getCellInfo(ed,b,cell_k)
            
            cellInfo = ed.cellInfo{cell_k};
            baseDir = ed.baseDirs{b};
            batNum = ed.batNums{b};
            
            
            
            switch ed.expType
                
                case 'adult'
                    
                    spikeDir = fullfile(baseDir,'spike_data',[batNum '_' cellInfo '.csv']);
                    %       If we want to calculate spike waveform stats, need to include path to actual .ntt file
                    ttDir = [];
                    call_data_dir = fullfile(baseDir,'call_data');
                    exp_date_str = cellInfo(1:8);
                    s = load(fullfile(ed.analysisDir,'cell_stability_info.mat'));
                    cell_stability_info = s.cell_stability_info;
                    idx = find(strcmp(cellInfo,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
                    stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
                    
                    if strcmp(ed.call_echo,'call') % are we looking at calls or echolocation clicks?
                        s = load(fullfile(call_data_dir,[exp_date_str '_cut_call_data.mat']));
                        cut_call_data = s.cut_call_data;
                        batIdx = strcmp({cut_call_data.batNum},batNum);
                        cut_call_data = cut_call_data(batIdx);
                        
                    elseif strcmp(ed.call_echo,'echo')
                        error('echolocation not yet implemented for adult data')
                    end
                    
                    audio2nlg = load(fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit.mat'])); % load fit data to sync audio to nlg data
                    stabilityBounds = stabilityBounds - audio2nlg.first_nlg_pulse_time;
                    success = true;
                    
                case 'adult_wujie'
                    
                    expDate = cellInfo(1:strfind(cellInfo,ed.tetrodeStr)-1);
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
                        batIdx = unique(cellfun(@(call) find(cellfun(@(bNum) strcmp(bNum,ed.batNums{b}),call)),{cut_call_data.batNum}));
                        
                        if length(batIdx) == 1
                            callpos = horzcat(cut_call_data.corrected_callpos);
                            callpos = callpos(batIdx,:);
                            [cut_call_data.corrected_callpos] = deal(callpos{:});
                        else
                            keyboard
                        end
                        call_info_fname = dir([audioDir 'call_info_*_' ed.call_echo '_' expDate '.mat']);
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
                    
                case 'juvenile'
                    
                    cluster_num_idx = strfind(cellInfo,ed.clusterStr)+length(ed.clusterStr);
                    tt_num_idx = strfind(cellInfo,ed.tetrodeStr )+length(ed.tetrodeStr );
                    
                    expDate = cellInfo(1:strfind(cellInfo,ed.tetrodeStr)-1);
                    tetrodeNum = str2double(cellInfo(tt_num_idx));
                    cellNum = str2double(cellInfo(cluster_num_idx:end));
                    
                    if cellNum < 10
                        sorted_cell_string = [ed.clusterStr '0'];
                    else
                        sorted_cell_string = ed.clusterStr;
                    end
                    
                    tetrode = [expDate ed.tetrodeStr num2str(tetrodeNum) sorted_cell_string num2str(cellNum)]; % build tetrode/cell string
                    ttDir = [ed.spike_data_dir 'bat' batNum filesep expDate filesep tetrode '.ntt']; % filename for cell of interest
                    audioDir = [baseDir 'bat' batNum filesep 'neurologger_recording' expDate '\audio\ch1\']; % directory where .wav files are stored (and has a subfolder 'Analyzed_auto')
                    s = load([ed.analysisDir 'cell_stability_info']);
                    cell_stability_info = s.cell_stability_info;
                    idx = find(strcmp(tetrode,{cell_stability_info.cellInfo}) & strcmp(batNum,{cell_stability_info.batNum}));
                    stabilityBounds = [cell_stability_info(idx).tsStart cell_stability_info(idx).tsEnd];
                    
                    
                    try
                        
                        if strcmp(ed.call_echo,'call') % are we looking at calls or echolocation clicks?
                            s = load([audioDir 'cut_call_data.mat']);
                            cut_call_data = s.cut_call_data;
                            cut_call_data = cut_call_data(~[cut_call_data.noise]);
                            
                        elseif strcmp(ed.call_echo,'echo')
                            s = load([audioDir 'cut_echo_data.mat']);
                            cut_call_data = s.cut_call_data;
                            cut_call_data = cut_call_data(~[cut_call_data.noise]);
                            
                            s = load([audioDir 'juv_call_info_' ed.batNums{b} '_echo.mat']);
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
                    spikeDir = [];
            end
            
        end
        
    end
    
end