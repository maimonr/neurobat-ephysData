classdef ephysData
    properties
        baseDirs
        batNums
        boxNums
        dateFormat
        birthDates
        analysisDir
        activeChannels
        spike_data_dir
        lfp_data_dir
        lfp_fs
        expType
        callType
        
        serverStr = 'server1';
        pathStr = 'users\maimon\';
        
        remote_drive_letter
        serverPath
        recLogs
    end
    
    methods
        
        function ed = ephysData(expType)
            
            if nargin == 0
                ed.expType = 'juvenile';
            else
                ed.expType = expType;
            end
            
            ed.remote_drive_letter = get_remote_drive_letter(ed);
            ed.serverPath = get_server_path(ed);
            
            switch ed.expType
                
                case 'juvenile'
                    
                    ed.batNums = {'71319','71284','71173','13681'};
                    ed.baseDirs = repmat({'E:\ephys\juvenile_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = {datetime(2016,4,23),datetime(2016,09,24),datetime(2016,09,21),datetime(2017,10,20)};num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = repmat({'C:\Users\phyllo\Documents\Maimon\ephys\data_analysis_results\'},1,length(ed.batNums));
                    ed.activeChannels = {setdiff(0:15,[1 4]),setdiff(0:15,[1:4, 8, 9, 15]),setdiff(0:15,[12 13]),setdiff(0:15,8)};
                    ed.spike_data_dir = repmat({'E:\ephys\juvenile_recording\tetrode_data\'},1,length(ed.batNums));
                    ed.lfp_fs = 29297;
                    ed.lfp_data_dir = [];
                    
                case 'adult'
                    
                    ed.batNums = {'14620','71334','65694','71360'};
                    ed.boxNums = num2cell(nan(1,length(ed.batNums)));
                    ed.baseDirs = repmat({'E:\ephys\adult_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = repmat({'E:\ephys\adult_recording\data_analysis_results\'},1,length(ed.batNums));
                    ed.activeChannels = {0:15,0:15,0:15,[4:7 12:15]};
                    ed.spike_data_dir = repmat({'E:\ephys\adult_recording\spike_data\'},1,length(ed.batNums));
                    ed.lfp_fs = 31250;
                    ed.lfp_data_dir = [];
                    
                case 'adult_operant'
                    
                    ed.batNums = {'59886','65705','71382','65702'};
                    ed.boxNums = {'1','1','2','2'};
                    ed.baseDirs = repmat({'E:\ephys\adult_operant_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = repmat({'E:\ephys\adult_operant_recording\data_analysis_results\'},1,length(ed.batNums));
                    ed.activeChannels = {setdiff(0:15,10),0:15,0:15,setdiff(0:15,5)};
                    ed.spike_data_dir =  repmat({'E:\ephys\adult_operant_recording\spike_data\'},1,length(ed.batNums));
                    ed.lfp_fs = 31250;
                    ed.lfp_data_dir = [];
                    
                case 'adult_social'
                    
                    ed.batNums = {'13688','11636','14612','11682','71216'};
                    ed.boxNums = num2cell(nan(1,length(ed.batNums)));
                    
                    ed.baseDirs = repmat({ed.serverPath},1,length(ed.batNums));
                    
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = repmat({fullfile(ed.serverPath,'data_analysis_results')},1,length(ed.batNums));
                    ed.activeChannels = {0:15,0:15,setdiff(0:15,8:9),0:15,0:15};
                    ed.spike_data_dir =  repmat({fullfile(ed.serverPath,'spike_data')},1,length(ed.batNums));
                    ed.lfp_fs = 31250;
                    ed.lfp_data_dir = repmat({fullfile(ed.serverPath,'lfp_data')},1,length(ed.batNums));
                    
                case 'adult_wujie'
                    
                    ed.batNums = {'59813','71348','65997'};
                    ed.baseDirs = repmat({'E:\ephys\adult_recording\wujie_data\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = 'E:\ephys\wujie_adult_recording\data_analysis_results\';
                    ed.activeChannels = [];
                    ed.spike_data_dir = 'E:\ephys\wujie_adult_recording\spike_data\';
                    ed.lfp_data_dir = 'E:\ephys\wujie_adult_recording\lfp_data\';
                    
            end
            ed.recLogs = readtable(fullfile(ed.baseDirs{1},'documents','recording_logs.csv'));
            
        end
    end
end

function remote_drive_letter = get_remote_drive_letter(ed)
[~,caption_str]= dos('wmic logicaldisk get caption');
[~,name_str]= dos('wmic logicaldisk get volumename');

caption_str = strsplit(caption_str,'\n');
name_str = strsplit(name_str,'\n');

remote_drive_letter = caption_str{contains(name_str,ed.serverStr)};
remote_drive_letter = deblank(remote_drive_letter);
end

function serverPath = get_server_path(ed)
serverPath = fullfile(ed.remote_drive_letter,ed.pathStr,[ed.expType '_recording']);
end