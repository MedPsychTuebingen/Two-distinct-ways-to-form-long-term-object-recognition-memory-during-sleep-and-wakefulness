%% SO phase - spindle amplitude coupling for frontal EEG elecs 
% Phase amplitude coupling for the duration of SO

clear
%% Paths
dirProject      = 'Y:\Max\1wk_NOR_new\';
addpath(strcat(dirProject, 'Scripts'));
addpath('Y:\Max\CircStat2012a');

%% Basics
fsample         = 1000;
recLength       = 120; %in min

% List of EEG Channels
numChannels     = 2;
strChannel      = cell(1,numChannels);
strChannel{1,1} = 'EEG_Left';
strChannel{1,2} = 'EEG_Right';

% Filters
SO_LowPassFilter   = 4.5; 
SO_HighPassFilter  = 0.3; 
Spi_LowPassFilter  = 16;
Spi_HighPassFilter = 10;

% EEG filter SO band (narrow for phase extraction)
LowPassFilter_SO_Phase  = 1.25;
HighPassFilter_SO_Phase = 0.5;
        
%% Loading RecordingInfo & Preparing Variables
Info = readtable(strcat(dirProject,'Info_1wk_NOR.xlsx'));

%% Coocurrence analysis
for iAnimal = 1:size(Info,1)
    for iCh = 1:numChannels

        % Percentage of co-occurance of detected SO and spindles 
        % Read in Detections
        DetectedSO  = load(string(strcat('Detections\SO_', Info.Name(iAnimal, 1), '.mat')), strChannel{1,iCh});
        DetectedSO  = DetectedSO.(strcat(strChannel{1,iCh})).trialinfo;
        DetectedSpi = load(string(strcat('Detections\Spi_', Info.Name(iAnimal, 1), '.mat')), strChannel{1,iCh});
        DetectedSpi = DetectedSpi.(strcat(strChannel{1,iCh})).trialinfo;
        
        % Read EEG data
        tmpEEG = load(strcat(dirProject,'MATLABFiles\',Info.Name{iAnimal,1},'.mat'), strcat('Ch',num2str(Info{iAnimal,4+iCh})));
        tmpEEG = tmpEEG.(strcat('Ch',num2str(Info{iAnimal,4+iCh})));
        tmpEEG.values = tmpEEG.values(1:recLength*60*fsample,:);
        tmpEEG.length = recLength*60*fsample;
        tmpEEG.values = tmpEEG.values*-1;  % Invert the data    
        
        % Location of downstate peak, Begin and End of SOs 
        SONegPkLoc = DetectedSO(:,5); % can be read directly from detections
        SOBegin    = DetectedSO(:,2); 
        SOEnd      = DetectedSO(:,3);
       
        % Location of maximum spindle amplitude
        % Spindle band filters 
        [d,e]           = butter(3,2*Spi_LowPassFilter/fsample,'low'); 
        FilteredEEG_Spi = filtfilt(d,e,tmpEEG.values); 
        [d,e]           = butter(3,2*Spi_HighPassFilter/fsample,'high'); 
        FilteredEEG_Spi = filtfilt(d,e,FilteredEEG_Spi);
        SpiEnvelope     = abs(hilbert(FilteredEEG_Spi)); 
        
        SpiMaxAmp = zeros(length(DetectedSpi(:,1)),2);
        j = 1;
        for iSpi = 1:size(DetectedSpi)
            [MaxAmp, Idx] = max(SpiEnvelope(DetectedSpi(iSpi,2):DetectedSpi(iSpi,3)));
            SpiMaxAmp(j+length(MaxAmp)-1,1) = DetectedSpi(iSpi,2) + Idx;
            SpiMaxAmp(j+length(MaxAmp)-1,2) = MaxAmp;
            j = j+length(MaxAmp);
        end
        clear j Idx MaxAmp iSpi;
        
        % Narrow SO band filters
        [d,e]          = butter(3,2*LowPassFilter_SO_Phase/fsample,'low');
        FilteredEEG_SO_Phase    = filtfilt(d,e,tmpEEG.values); 
        [d,e]                   = butter(3,2*HighPassFilter_SO_Phase/fsample,'high'); 
        FilteredEEG_SO_Phase    = filtfilt(d,e,FilteredEEG_SO_Phase);
        EEGHilbert_SO_Phase     = hilbert(FilteredEEG_SO_Phase); % Hilbert 
        SO_phase                = angle(EEGHilbert_SO_Phase); % extract phase in radians
        
        % Calculate co-occurrance (+/- 1.2 s window) 
        SO_Spi = zeros(length(DetectedSO),5);
        j = 1;
        for iSO = 1:length(SONegPkLoc) 
            % Find Phase = 90 degrees in falling flank (start of SO)
            falling_flank = rad2deg(SO_phase(SONegPkLoc(iSO)-1000: SONegPkLoc(iSO)));
            [~,idx]=min(abs(90-falling_flank));
            begin_idx = SONegPkLoc(iSO)-1000 + idx;
            clear falling_flank idx;
            
            % Find Phase = 90 degrees in rising flank (end of SO)
            if SONegPkLoc(iSO)+1500 > 7200000 % If SO at end of recording set end of SO phase to end of recording
                end_idx = 7200000;
            else
                rising_flank = rad2deg(SO_phase(SONegPkLoc(iSO): SONegPkLoc(iSO)+1500));
                [~,idx]=min(abs(90-rising_flank));
                end_idx = SONegPkLoc(iSO) + idx;
                clear idx;
            end
           
            % Find Spindles that are nested in detected SO (within complete
            % SO cycle)
            tmpSO_Spi = find(begin_idx <= SpiMaxAmp & SpiMaxAmp <= end_idx);
                    
            if ~isempty(tmpSO_Spi)
                SO_Spi(j:j+length(tmpSO_Spi)-1,1) = tmpSO_Spi(:,1); % stores position in SpiMinLoc of spindle cooccurring with SO 
                SO_Spi(j:j+length(tmpSO_Spi)-1,2) = SONegPkLoc(iSO); % stores location of downstate of SO 
                SO_Spi(j:j+length(tmpSO_Spi)-1,3) = DetectedSO(iSO,2); % stores begin of SO 
                SO_Spi(j:j+length(tmpSO_Spi)-1,4) = DetectedSO(iSO,3); % stores end of SO 
                SO_Spi(j:j+length(tmpSO_Spi)-1,5) = iSO; % SO number
                j = j+length(tmpSO_Spi);
            else
                clear tmpSO_Spi;
            end
        end
        
        SO_Spi = SO_Spi(any(SO_Spi ~= 0,2),:); % removes all zeros
        SO_Spi(:,6) = DetectedSpi(SO_Spi(:,1),2); % stores begin of spindle
        SO_Spi(:,7) = DetectedSpi(SO_Spi(:,1),3); % stores end of spindle  
        SO_Spi(:,8) = SpiMaxAmp(SO_Spi(:,1),1); % position of max Amp of Spindle 

        %% Phase (SO) - Amplitude (spindle) coupling 
        
        % Extract SO phase for SO-Spi coupled events
        Spi_excluded = 0;
        for iSO = 1:size(SO_Spi,1)
            % Find Phase = 90 degrees in falling flank (start of SO)
            falling_flank = rad2deg(SO_phase(SO_Spi(iSO, 2)-1000: SO_Spi(iSO, 2)));
            [~,idx]=min(abs(90-falling_flank));
            begin_idx = SO_Spi(iSO, 2)-1000 + idx;
            clear falling_flank idx;
            
            % Find Phase = 90 degrees in rising flank (end of SO)
            rising_flank = rad2deg(SO_phase(SO_Spi(iSO, 2): SO_Spi(iSO, 2)+1500));
            [~,idx]=min(abs(90-rising_flank));
            end_idx = SO_Spi(iSO, 2) + idx;
            clear idx;
            
            if SO_Spi(iSO,8)-begin_idx < 0 || SO_Spi(iSO,8)-end_idx > 0
                Spi_excluded = Spi_excluded + 1; % exclude Spindles from analysis if Max Amp is outside of begin/end of SO phase
                SO_Spi(iSO,9) = nan;
            else
                tmp_SO_phase = SO_phase(begin_idx:end_idx);
                SO_Spi(iSO,9) = tmp_SO_phase(SO_Spi(iSO,8)-begin_idx); % store SO phase at Max Spindle Amp
            end
        end
        
        %% Save output 
        % Event info 
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh)))) = SO_Spi;
        
        % Number of spindle centers occurring within a +/- 1.2 s window
        % around the downstate peak of a SO, expressed as the ratio of 
        % Co-occuring events 
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,1) =  size(SO_Spi,1); 
        
        % Co-occuring SOs
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,2) =  length(unique(SO_Spi(:,5))); 
        
        % Percent of coupled SOs
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,3) = length(unique(SO_Spi(:,5)))/size(DetectedSO,1)*100;  
        
        % Co-ocurring Spindles
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,4) = length(unique(SO_Spi(:,1)));  
        
        % Percent of coupled SOs
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,5) = length(unique(SO_Spi(:,1)))/size(DetectedSpi,1)*100;
        
        % preferred SO angle for spindle (mean) 
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,6) = circ_mean(SO_Spi(:,9));  
        
        % preferred SO angle for spindle (median)
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,7) = circ_median(SO_Spi(:,9));  
        
        % Number of excluded spindles for Phase-Amplitude Coupling Analysis
        Coupling.(strcat('SO_Spi_',Info.Name{iAnimal,1}(1:2),Info.Name{iAnimal,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh), '_sum')))(1,8) = Spi_excluded;      
        
        %clear SO_Spi SpiMaxAmp tmpEEG SpiEnvelope SONegPkLoc SO_phase Spi_excluded
    end
end

%% Save data
for i = 1:length(Info.Name)
    file = Coupling.(strcat('SO_Spi_',Info.Name{i,1}(1:2),Info.Name{i,1}(end-1:end)));
    save(strcat(dirProject,'Coupling_narrow\SO_Spi_',Info.Name{i,1}(1:2),Info.Name{i,1}(end-1:end),'.mat'), '-struct', 'file' ,'-v6')
    clear file
end

%% Plot single histos for all animals 
for iRec = 1:size(Info,1)
    for iCh = 1:numChannels
    subplot(4,6,iRec)
    tmp_angle =  Coupling.(strcat('SO_Spi_',Info.Name{iRec,1}(1:2),Info.Name{iRec,1}(end-1:end)))...
            .(string(strcat(strChannel(1,iCh))))(:,9);
    circ_plot(tmp_angle, 'hist', []);
    title(string(strcat(Info.Name{iRec,1}(1:2),Info.Name{iRec,1}(end-1:end))));
    hold all
    clear tmp_angle;
    end
end

close all;
