classdef ephysData
    properties (Constant = true)
        baseDirs = {'Z:\users\Maimon\ephys\','E:\ephys\juvenile_recording\','E:\ephys\juvenile_recording\'};
        batNums = {'71319','71284','71173'};
        dateFormat = 'yyyyMMdd';
        birthDates = {datetime(2016,4,23),datetime(2016,09,24),datetime(2016,09,21)};
        analysisDir = 'C:\Users\phyllo\Documents\Maimon\ephys\data_analysis_results\';
    end
    properties
        call_echo = 'call'
        
    end
    
end