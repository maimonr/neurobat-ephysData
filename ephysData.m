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
        expType
        callType = 'call'
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
                    ed.activeChannels = {0:15,0:15,0:15,[5:8 13:16]};
                    ed.spike_data_dir = 'E:\ephys\adult_recording\spike_data\';
                    ed.lfp_data_dir = [];
                    
                case 'adult_operant'
                    
                    ed.batNums = {'59886','65705','71382','65702'};
                    ed.boxNums = {'1','1','2','2'};
                    ed.baseDirs = repmat({'E:\ephys\adult_operant_recording\'},1,length(ed.batNums));
                    ed.dateFormat = 'yyyyMMdd';
                    ed.birthDates = num2cell(repmat(NaT,1,length(ed.batNums)));
                    ed.analysisDir = 'E:\ephys\adult_operant_recording\data_analysis_results\';
                    ed.activeChannels = {setdiff(0:15,10),0:15,0:15,setdiff(0:15,5)};
                    ed.spike_data_dir = 'E:\ephys\adult_operant_recording\spike_data\';
                    ed.lfp_data_dir = [];
                    
                    
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
            
        end
        
    end
    
end