classdef ephysData
    properties (Constant = true)
        baseDirs = {'E:\ephys\juvenile_recording\','E:\ephys\juvenile_recording\','E:\ephys\juvenile_recording\','E:\ephys\juvenile_recording\'};
        batNums = {'71319','71284','71173','13681'};
        dateFormat = 'yyyyMMdd';
        birthDates = {datetime(2016,4,23),datetime(2016,09,24),datetime(2016,09,21),datetime(2017,10,20)};
        analysisDir = 'C:\Users\phyllo\Documents\Maimon\ephys\data_analysis_results\';
        activeChannels = {setdiff(0:15,4),setdiff(0:15,[1:4, 8, 9]),setdiff(0:15,12),setdiff(0:15,8)};
    end
    properties
        call_echo = 'call'
        
    end
    
end