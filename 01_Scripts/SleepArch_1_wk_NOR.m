%% Basic sleep parameters from sleep scoring

clear
%% Paths
dirProject      = 'Y:\Max\1wk_NOR_new\';
addpath(strcat(dirProject, 'Scripts'));

%% Basics
fsample         = 1000;
recLength       = 120; %in min

%% Loading RecordingInfo
Info = readtable(strcat(dirProject,'Info_1wk_NOR.xlsx'));

for iRecs = 1:size(Info,1)
    
    tmpRawSleepScoring  = dlmread(strcat(dirProject, 'SleepScoring/', Info.Name{iRecs,1},'_scoring.txt'),'\t',[1 0 0 2]);
    
    %%% Checking length
    if size(tmpRawSleepScoring,1)*10/60 >= 120
        tmpRawSleepScoring = tmpRawSleepScoring(1:recLength*60/10,:);
        %         %%% Analysis 120 min
        %         SleepArch.(strcat('A_',Info.Name{iRecs,1}(1:2),Info.Name{iRecs,1}(end-1:end))) = SleepArchFunc(tmpRawSleepScoring,recLength*60*fsample,fsample);
        %%% Analysis 120 min
        for iTimeBin = 1
            if iTimeBin == 1
                tmpBeg = 1;
                tmpEnd = recLength*60/10;
            end
            tmpRawSleepScoring_bin = tmpRawSleepScoring(tmpBeg:tmpEnd,:);
            SleepArch.(strcat('Arch_',Info.Name{iRecs,1}(1:2),Info.Name{iRecs,1}(end-1:end))){1,iTimeBin} = SleepArchFunc(tmpRawSleepScoring_bin,recLength*60*fsample,fsample);
        end
        
%         %%% Analysis per Hour
%         for iTimeBin = 1:2
%             if iTimeBin == 1
%                 tmpBeg = 1;
%                 tmpEnd = 60*60/10;
%             else
%                 tmpBeg = 60*60/10+1;
%                 tmpEnd = recLength*60/10;
%             end
%             tmpRawSleepScoring_bin = tmpRawSleepScoring(tmpBeg:tmpEnd,:);
%             SleepArch.(strcat('A_',Info.Name{iRecs,1})){1,iTimeBin} = SleepArchFunc(tmpRawSleepScoring_bin,recLength*60*fsample,fsample);
%         end
    else
        error('Error. Recording shorter than 2 hours')
    end
    
    clear tmpRawSleepScoring
end

% save(strcat(dirProject,'SleepArch\1_wk_NOR_SleepArch_120.mat'), '-struct', 'SleepArch','-v7.3')
save(strcat(dirProject,'SleepArch\1_wk_NOR 2020_SleepArch_120min.mat'), '-struct', 'SleepArch','-v7.3')