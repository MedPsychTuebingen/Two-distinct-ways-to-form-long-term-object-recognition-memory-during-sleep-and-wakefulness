%% Theta analysis for REM sleep epochs for frontal EEG electrodes 
clear
%% Paths
dirProject      = 'add_path';
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

% Filters for theta freq band 
LowPassFilter   = 10; 
HighPassFilter  = 5; 

%% Extract Theta Power for each REM epoch 
for iRecs = 1%:size(Info,1)
    tmpRawSleepScoring  = dlmread(strcat(dirProject, 'SleepScoring/', Info.Name{iRecs,1},'_scoring.txt'),'\t',[1 0 0 2]);
    
    %%% Cutting the Scoring to 120 min
    tmpRawSleepScoring  = tmpRawSleepScoring(1:recLength*60/10,:);
    
    % Deleting Epochs with Artefacts
    tmpRawSleepScoring((tmpRawSleepScoring(:,3)==1),2) = 0;
    
    % Extend RawScoring to the length of the recording
    tmpSleepScoring     = zeros((length(tmpRawSleepScoring)*fsample*10),1);
    
    % Extending the recording
    for iExt = 1:(length(tmpRawSleepScoring)*fsample*10)
        tmpSleepScoring(iExt,1)     = tmpRawSleepScoring(ceil((1/fsample)*iExt/10),2);
    end
    clear iExt tmpRawSleepScoring
    
    % Find REM epochs
    REM = find (tmpSleepScoring(:,1)==3); %find REM episodes
    
    if isempty(REM) % in case there is no REM sleep in entire recording
        REMEpisodes = [];
    else
        REMEndEpisode = [];
        REMBegEpisode = [];
        
        for i=2:length(REM)-1
            if REM(i) - REM(i-1) > 1
                REMBegEpisode = [REMBegEpisode,REM(i)];
            end
            if REM(i+1) - REM(i) > 1
                REMEndEpisode = [REMEndEpisode,REM(i)];
                
            end
        end
        
        REMBegEpisode   = [REM(1),REMBegEpisode];
        REMEndEpisode   = [REMEndEpisode, REM(end)];
        
        REMEpisodes     = [REMBegEpisode;REMEndEpisode];
        %REMEpisodesBins.30mins = ;
        clear REMBegEpisode REMEndEpisode
    end
      
    % Filter signal
    for iCh = 1:numChannels %size(whos('-file',strcat(dirProject,'MATLABFiles/', Info.Name{iRecs,1},'.mat')),1) 
        tmpEEG = load(strcat(dirProject,'MATLABFiles/',Info.Name{iRecs,1},'.mat'), strcat('Ch',num2str(Info{iRecs,4+iCh})));
        tmpEEG = tmpEEG.(strcat('Ch',num2str(Info{iRecs,4+iCh})));
        Channel = tmpEEG.title;
        % Invert the data
        tmpEEG.values = tmpEEG.values*-1;     

        [d,e]               = butter(3,2*LowPassFilter/fsample,'low'); %Use butterworth filter 3rd order. 
        FiltSignal          = filtfilt(d,e,tmpEEG.values); %Filter Signal
        [d,e]               = butter(3,2*HighPassFilter/fsample,'high'); %Use butterworth filter 3rd order. 
        FiltSignal          = filtfilt(d,e,FiltSignal);
        FiltEnv             = abs(hilbert(FiltSignal));
        
        %data.(strcat('Ch',num2str(Info{iRecs,4+iCh}))).FiltEnv = FiltEnv;
        clear tmpEEG FiltSignal d e
        
        % Calculate Power for each REM epoch
         tmpREMtrialinfo.(strcat(Channel)) = zeros(length(REMEpisodes),3);
         tmpREMtrialinfo.(strcat(Channel, '_Hilbert')) = zeros(7200000,1);
         tmpREMtrialinfo.(strcat(Channel, '_sum')) = zeros(1,4);
         j = 1;
         for iREMTrial = 1:size(REMEpisodes,2) 
             
            tmpREMEnv       = FiltEnv(REMEpisodes(1,iREMTrial):REMEpisodes(2,iREMTrial));
            tmpREMEnergy     = sum(tmpREMEnv);
            
            % Organizing Trialinfo
            tmpREMtrialinfo.(strcat(Channel))(iREMTrial,1) = (ceil(REMEpisodes(1,iREMTrial)/1800000)*30); % time bin in 30 mins
            tmpREMtrialinfo.(strcat(Channel))(iREMTrial,2) = (diff(REMEpisodes(:,iREMTrial))+1)/fsample; %Duration in Secs
            tmpREMtrialinfo.(strcat(Channel))(iREMTrial,3) = tmpREMEnergy; % Energy 
            tmpREMtrialinfo.(strcat(Channel, '_Hilbert'))(j:j+length(tmpREMEnv)-1,1) = tmpREMEnv; % absolute Hilbert of REM epochs
            j = j+length(tmpREMEnv);
            clear tmpREMPower
         end
         tmpREMtrialinfo.(strcat(Channel, '_Hilbert')) = tmpREMtrialinfo.(strcat(Channel, '_Hilbert'))...
             (any(tmpREMtrialinfo.(strcat(Channel, '_Hilbert')) ~= 0,2),:); % removes all zeros
         
         % Calculate means
         tmpREMtrialinfo.(strcat(Channel, '_sum'))(1,1) = length(tmpREMtrialinfo.(strcat(Channel, '_Hilbert')))/fsample/60; % Total duration REM in mins
         tmpREMtrialinfo.(strcat(Channel, '_sum'))(1,2) = mean(tmpREMtrialinfo.(strcat(Channel))(:,2)); % Mean duration REM epochs in seconds
         tmpREMtrialinfo.(strcat(Channel, '_sum'))(1,3) = mean(tmpREMtrialinfo.(strcat(Channel, '_Hilbert'))); % Mean power 
         tmpREMtrialinfo.(strcat(Channel, '_sum'))(1,4) = sum(tmpREMtrialinfo.(strcat(Channel))(:,3)); % Total energy 
         tmpREMtrialinfo = rmfield(tmpREMtrialinfo, (strcat(Channel, '_Hilbert')));
    end    
 
    clear tmpREMEpisode
    
    REMPower.(strcat('REMTheta_',Info.Name{iRecs,1}(1:2),Info.Name{iRecs,1}(end-1:end))) = tmpREMtrialinfo;
    %clear tmpREMtrialinfo
end

% Save data
for i = 1:length(Info.Name)
    file = REMPower.(strcat('REMTheta_', Info.Name{i,1}(1:2),Info.Name{i,1}(end-1:end)));
    save(strcat(dirProject,'REMTheta\REMTheta_', Info.Name{i,1}(1:2),Info.Name{i,1}(end-1:end),'.mat'), '-struct', 'file' ,'-v6')
    clear file
end

