%% Spindle and SO detections for frontal EEG screw electrodes
clear
%% Paths
dirProject      = 'Y:\Max\1wk_NOR_new\';
addpath(strcat(dirProject, 'Scripts'));
%% Basics
fsample         = 1000;

% List of EEG Channels
numChannels     = 2;
strChannel      = cell(1,numChannels);
strChannel{1,1} = 'EEG_Left';
strChannel{1,2} = 'EEG_Right';

%% Analysis
recLength = 120; %in min

%% Loading RecordingInfo & Preparing Variables
Info = readtable(strcat(dirProject,'Info_1wk_NOR.xlsx'));

for iRecs = 1:size(Info,1)
    tmpRawSleepScoring  = dlmread(strcat(dirProject, 'SleepScoring/', Info.Name{iRecs,1},'_scoring.txt'),'\t',[1 0 0 2]);
    
    %%% Cutting the Scoring to 120 min
    tmpRawSleepScoring  = tmpRawSleepScoring(1:recLength*60/10,:);
    
    % Extending RawScoring to the length of the recording
    tmpSleepScoring     = zeros((length(tmpRawSleepScoring)*fsample*10),1);
    tmpSleepArtefacts   = zeros((length(tmpRawSleepScoring)*fsample*10),1);
    
    % Extending the recording
    for iExt = 1:(length(tmpRawSleepScoring)*fsample*10)
        tmpSleepScoring(iExt,1)     = tmpRawSleepScoring(ceil((1/fsample)*iExt/10),2);
        tmpSleepArtefacts(iExt,1)   = tmpRawSleepScoring(ceil((1/fsample)*iExt/10),3);
    end
    clear iExt tmpRawSleepScoring
    
    for iCh = 1:numChannels%size(whos('-file',strcat('MATLABFiles/', Info.Name{iRecs,1},'.mat')),1) %-1, because of EMG
        tmpEEG = load(strcat(dirProject,'MATLABFiles/',Info.Name{iRecs,1},'.mat'), strcat('Ch',num2str(Info{iRecs,4+iCh})));
        tmpEEG = tmpEEG.(strcat('Ch',num2str(Info{iRecs,4+iCh})));
        
        % Invert the data
        tmpEEG.values = tmpEEG.values*-1;      
        
        %% Event Detection
        tmpspidetection = SpindleDetection_1wk_NOR(tmpEEG.values,tmpSleepScoring,tmpSleepArtefacts,tmpEEG.title,fsample);
        tmpSOdetection  = SODetection_1wk_NOR(tmpEEG.values,tmpSleepScoring,tmpSleepArtefacts,tmpEEG.title,fsample);
        
        % Storing output
        SODetection.(tmpEEG.title)          = tmpSOdetection;
        SpindleDetection.(tmpEEG.title)     = tmpspidetection;
        clear tmpSOdetection tmpspidetection tmpEEG
    end
    clear tmpSleepScoring tmpSleepArtefacts
    
    %% Saving output per Animal
    save(strcat(dirProject,'Detections/SO_',Info.Name{iRecs,1},'.mat'),'-struct', 'SODetection','-v7.3')
    save(strcat(dirProject,'Detections/Spi_',Info.Name{iRecs,1},'.mat'),'-struct', 'SpindleDetection','-v7.3')
end