function out = SpindleDetection_1wk_NOR(RecordingValues, SleepScoring, ScoringArtefacts, Label, fsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindle Detection Algorith for EEG/LFP siganls
% INPUT:    RecordingValues: Matrix containing Raw EEG or LFP values
%           SleepScoring: Matrix containing Sleep Scoring for each data point
%           ScoringArtefact: Matrix coding the artefacts (1) of the signal
%           for each data point.
%           Label: CHAR with the name of the Channel used.
%           fsample: sampling rate.
% OUT:      Output contains detection info in the following way:
%           out.label           = CELL containg the name of the Channel(s)
%           out.fsample         = sampling rate
%           out.trial           = Event(s) raw signal
%           out.time            = time/length of the Trial(s)
%           out.trialinfo       = Detailed info about each one of the trials.
%               out.trialinfo(:,1) = To which time bin corresponds that particular event 
%               out.trialinfo(:,2) = startTime: begin of event 
%               out.trialinfo(:,3) = endTime: spindle: negative zero crossing
%               out.trialinfo(:,4) = duration: duration from start to end in seconds (spindle: betweet the two threshild crossings)
%               out.trialinfo(:,5) = minTime: time of minimum (spindle: largest negative trough during spindle) in datapoints of original dataset
%               out.trialinfo(:,6) = minAmp: amplitude of minimum (spindle: largest negative trough during spindle) in �V
%               out.trialinfo(:,7) = Power
%               out.trialinfo(:,8) = Frequency
%           out.vector          = Vector containing the presence of the events
%           out.sleeptimeinsec  = Time in SECONDS used for the Detections (excl. time spent in Artefacts)
%           out.density_permin  = Density of events for the whole recording time (excl. time spent in Artefacts)
%
% Authors:  Carlos N. Oyanedel - Uni T�bingen, Germany
%           Niels Niethard - Uni T�bingen, Germany
%           Thanks to Dr. Hong-Viet Ngo, University of Birmingham, UK
% Contact:  jan.born@uni-tuebingen.de
%
% Date:     21.08.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic Parameters
% Spindle Duration in Seconds
MinSpindleDur = 0.4;
MaxSpindleDur = 2;

% Filters in Hz
LowPassFilter       = 16;
HighPassFilter      = 10;

% Sleep Stage of Interest
%1: WAKE 2: NREM; 3:REM; 4:PREM 
SleepStages     = 1:4;
SleepStageToUse = 2;

% STD above mean
StdValue        = 1.5; %The amount of STD above the mean for the detection

% Trial Window in seconds
TrialWindow  = 2.5; % Window before and after the Min Peak

%% Filtering the signal and calculating the envelope
[d,e]               = butter(3,2*LowPassFilter/fsample,'low'); %Use butterworth filter 3rd order. 15Hz cutoff
FilteredEEG         = filtfilt(d,e,RecordingValues); %Filter Signal
[d,e]               = butter(3,2*HighPassFilter/fsample,'high'); %Use butterworth filter 3rd order. 8Hz cutoff
FilteredEEG         = filtfilt(d,e,FilteredEEG);
FilteredEEGEnvelope = abs(hilbert(FilteredEEG));

%% Smoothing the FilteredEEGEnvelope
run SmoothFilteredEEGEnvelope.m

%% Mean and STD of the SmoothedFilteredEEGEnvelope of an specific Sleep Stage
ValidSleepScoring                       = SleepScoring;
ValidSleepScoring(ScoringArtefacts==1)  = 0; %Delete data with artefacts
SmoothedFilteredEEGEnvelopeMean         = mean(SmoothedFilteredEEGEnvelope(find(ValidSleepScoring==SleepStageToUse),1));
SmoothedFilteredEEGEnvelopeSTD          = std(SmoothedFilteredEEGEnvelope(find(ValidSleepScoring==SleepStageToUse)));

%% SpindleVector.
SpindleVector   = zeros(length(RecordingValues),1);

for i=1:length(RecordingValues) %Detecting Spindles crossing the power threshold
    if SmoothedFilteredEEGEnvelope(i) > (SmoothedFilteredEEGEnvelopeMean+StdValue*SmoothedFilteredEEGEnvelopeSTD)
        SpindleVector(i) = 1;
    end
end

%% If want to delete event detected in different sleep stages
% Delete detected values in the other Sleep Stages
SpindleVector(1,1)                      = 0;
SleepStagetoRemove                      = SleepStages(SleepStages~=SleepStageToUse);
SpindleVector(ValidSleepScoring==SleepStagetoRemove(1)) = 0; %
SpindleVector(ValidSleepScoring==SleepStagetoRemove(2)) = 0; %
SpindleVector(ValidSleepScoring==SleepStagetoRemove(3)) = 0; %
SpindleVector(ValidSleepScoring==8)     = 0; %In Artefacts 8
SpindleVector(ValidSleepScoring==0)     = 0; %In case of mistmatch in extending the Sleep Scoring ,i.e. DiffScoringRec>0
SpindleVector(length(ValidSleepScoring)+1:end) = 0; 

%% Detecting Event Beg and End
% Position of values detected
SpindleVectorLoc    = find(SpindleVector==1); 

% Beg and End of each event.
SpindleEnd =[];
SpindleBeg =[];

for i=2:length(SpindleVectorLoc)-1
    if SpindleVectorLoc(i) - SpindleVectorLoc(i-1) > 1
        SpindleBeg = [SpindleBeg,SpindleVectorLoc(i)];
    end
    if SpindleVectorLoc(i+1) - SpindleVectorLoc(i) > 1
        SpindleEnd = [SpindleEnd,SpindleVectorLoc(i)];
    end
end

SpindleBeg          = [SpindleVectorLoc(1), SpindleBeg];
SpindleEnd          = [SpindleEnd, SpindleVectorLoc(end)];
SpindlesDetected    = [SpindleBeg; SpindleEnd];

%% Duration Threshold
% Now look for spindles with a Min and Max length
SpindleDuration = diff(SpindlesDetected);

MinSpindleDur = MinSpindleDur*fsample;
MaxSpindleDur = MaxSpindleDur*fsample;

j = 1;
for i = 1:size(SpindlesDetected,2)
    if SpindleDuration(i) <= MaxSpindleDur && SpindleDuration(i) >= MinSpindleDur %Checking which ones fulfill the duration criteria
        ValidSpindles(1,j)  = SpindlesDetected(1,i);
        ValidSpindles(2,j)  = SpindlesDetected(2,i);
        j                   = j+1;
    end
end
ValidSpindleDurationSec = (diff(ValidSpindles))./fsample;

% Peaks
PeakSpindleNumber = zeros(1,size(ValidSpindles,2));
for i = 1:length(ValidSpindles)
    PeakSpindleNumber(:,i) = length(findpeaks(FilteredEEG(ValidSpindles(1,i):ValidSpindles(2,i))));
end
% Frequency
tmpFreq = PeakSpindleNumber./ValidSpindleDurationSec;

%% Valid Ripple Vector
ValidSpindleVector = zeros(1,size(FilteredEEG,1));
for i=1:size(ValidSpindles,2)
    ValidSpindleVector(1,ValidSpindles(1,i):ValidSpindles(2,i)) = 1;
end

%% Spindle Power
tmppower = zeros(size(ValidSpindles,2),1);
for iPow = 1:size(ValidSpindles,2)
    tmppower(iPow,1) = trapz(SmoothedFilteredEEGEnvelope(ValidSpindles(1,iPow):ValidSpindles(2,iPow)));
end

%% Detecting Local Min
% Local Minima
MinValidSpindleVal = zeros(1,size(ValidSpindles,2)); % Local Minima Value Matrx
MinValidSpindleLoc = zeros(1,size(ValidSpindles,2)); % Local Minima Position
for iLocMin = 1:size(ValidSpindles,2)
    [MinValidSpindleVal(1,iLocMin), MinValidSpindleLoc(1,iLocMin)] = min(FilteredEEG(ValidSpindles(1,iLocMin):ValidSpindles(2,iLocMin)));
end

% Positions related to the whole recording for the Min and Max for each Valid Spindle Detected
ValidSpindleMinLoc = zeros(1,size(ValidSpindles,2));
for iPosMinMax = 1:size(ValidSpindles,2)
    ValidSpindleMinLoc(iPosMinMax)  = (ValidSpindles(1,iPosMinMax) + MinValidSpindleLoc(iPosMinMax));
end

%% Grand Average
SpindleGrandAverage = zeros(size(ValidSpindles,2),2*TrialWindow*fsample);
out.trial           = {};
out.time            = {};
WindowGrandAverage  = TrialWindow*fsample; % Window before and after the Min Peak
Cut                 = 0; %In case the Trial Window for the last detected event is longer than the recording 
CutBeg              = 0; %In case the Trial Window for the first detected starts before the recording file

for i = 1:size(ValidSpindles,2)
    if ((ValidSpindles(1,i)+MinValidSpindleLoc(i))+WindowGrandAverage-1) > size(RecordingValues,1)
        Cut = 1;
    elseif ((ValidSpindles(1,i)+MinValidSpindleLoc(i))-WindowGrandAverage-1) < 0
        CutBeg = 1;        
    else
        SpindleGrandAverage(i,:)    = RecordingValues((ValidSpindles(1,i)+MinValidSpindleLoc(i))-WindowGrandAverage:(ValidSpindles(1,i)+MinValidSpindleLoc(i))+WindowGrandAverage-1);
        out.trial{1,i}              = SpindleGrandAverage(i,:)*1000;
        out.time{1,i}               = -TrialWindow:1/fsample:TrialWindow-1/fsample;
    end
end

%% Creating eventdata - Output
out.label           = {Label};
out.fsample         = fsample;
tmptrialinfo        = zeros(size(ValidSpindles,2),13);
tmptrialinfo(:,1)   = (ceil(ValidSpindles(1,:)/1800000)*30); %To which timebin corresponds
tmptrialinfo(:,2)   = ValidSpindles(1,:); %startTime: begin of event (spindle: positive threshold crossing)
tmptrialinfo(:,3)   = ValidSpindles(2,:); %endTime: (spindle: negative zero crossing) 
tmptrialinfo(:,4)   = ValidSpindleDurationSec(:); %duration: duration from start to end in seconds (spindle: betweet the two threshild crossings)
tmptrialinfo(:,5)   = ValidSpindleMinLoc(:); %minTime: time of minimum (spindle: largest negative trough during spindle) in datapoints of original dataset
tmptrialinfo(:,6)   = MinValidSpindleVal(:)*1000; %minAmp: amplitude of minimum (spindle: largest negative trough during spindle) in �V
tmptrialinfo(:,7)   = tmppower'; %Power of each event detected
tmptrialinfo(:,8)   = tmpFreq'; %Freq of each event detected

if Cut == 1
   tmptrialinfo = tmptrialinfo(1:size(out.time,2),:);
end

if CutBeg == 1
    out.trial   = out.trial(2:end);
    out.time    = out.time(2:end);
end

out.trialinfo       = tmptrialinfo;
out.vector          = ValidSpindleVector;
out.sleeptime_sec  = ((size(find(SleepScoring == SleepStageToUse),1)) - (size(find(SleepScoring == SleepStageToUse & ScoringArtefacts == 1),1)))/fsample;
out.density_permin  = size(ValidSpindles,2)/(out.sleeptime_sec/60);
end