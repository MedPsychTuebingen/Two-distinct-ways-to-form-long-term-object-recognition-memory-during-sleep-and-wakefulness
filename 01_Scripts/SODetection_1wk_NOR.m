function out = SODetection_1wk_NOR(RecordingValues, SleepScoring, ScoringArtefacts, label, fsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SO Detection Algorith for EEG/LFP siganls
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
%               out.trialinfo(:,2) = startTime: SO: first up-down zero crossing
%               out.trialinfo(:,3) = endTime: spindle: negative zero crossing
%               out.trialinfo(:,4) = duration: duration from start to end in seconds (SO: between the two down-to-up crossings)
%               out.trialinfo(:,5) = minTime: time of minimum (SO: of negative half-wave/down-state) in datapoints of original dataset
%               out.trialinfo(:,6) = minAmp: amplitude of maximum (SO: of positive half-wave/up-state) in �V
%               out.trialinfo(:,7) = Amplitude
%               out.trialinfo(:,8) = Frequency
%           out.vector          = Vector containing the presence of the events
%           out.sleeptime_sec   = Time in SECONDS used for the Detections (excl. time spent in Artefacts)
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
% SO Duration in Seconds
MinSODur    = 0.8;
MaxSODur    = 2; 

% Filters in Hz
LowPassFilter   = 4.5; 
HighPassFilter  = 0.3; 

% Sleep Stage of Interest
%1: WAKE 2: NREM; 3:REM; 4:PREM 
SleepStageToUse = 2;

% Trial Window in seconds
TrialWindow  = 2.5; % Window: 2.5 sec before and after the Min Peak

%% Inverting the values for LFP signals
tmplabel = char(label);
tmplabel = tmplabel(1:3);
if strcmp(tmplabel(1:3),'CA1') == 1 || strcmp(tmplabel(1:3),'LFP') == 1
    RecordingValues = RecordingValues*-1;
end
    
%% Creating a artefact-free sleep scoring
ValidSleepScoring                       = SleepScoring;
ValidSleepScoring(ScoringArtefacts==1)  = 0;
ValidSleepScoring(SleepScoring==8)      = 0;

%% Filtering the signal;
[d,e]           = butter(3,2*LowPassFilter/fsample,'low'); %Use butterworth filter 3rd order. 15Hz cutoff
FilteredEEG     = filtfilt(d,e,RecordingValues); %Filter Signal
[d,e]           = butter(3,2*HighPassFilter/fsample,'high'); %Use butterworth filter 3rd order. 8Hz cutoff
FilteredEEG     = filtfilt(d,e,FilteredEEG);

%% Finding x-intercept from positive to negative
y                                       = diff(sign(FilteredEEG));
y(ValidSleepScoring~=SleepStageToUse)   = 0; %Excludes Sleep Stages that are not NREM, replacing the negative values for Zero, i.e it won't be considered in the indx.
y(length(SleepScoring)+1:end)           = 0;
indx_down                               = find(y<0); % crossing from positive to negative
indx_up                                 = find(y>0); % crossing from negative to positive

%% Calculating the difference to check the duration threshold
% First Check for the distribution of the duration of the detected crossing
% to negative.
% Hongi's values: 0.4-2s between x-crossing
MinSODur    = MinSODur*fsample; %Parameters according to the values we detected (Histogram)
MaxSODur    = MaxSODur*fsample; %Parameters according to the values we detected (Histogram)
j           = 1;
for i = 2:length(indx_down)
    SODiff = indx_down(i)-indx_down(i-1);
    if SODiff <= MaxSODur && SODiff >= MinSODur
        SOTemp1(1,j) = indx_down(i-1);
        SOTemp1(2,j) = indx_down(i);
        j=j+1;
    end
end

SOTemp1(1,:) = SOTemp1(1,:)+1; %Adding one position to the first position of the detected SO, in order to start with the first negative number, and not with the last postive before crossing zero.

%% Finding SO Peaks
% Local Maxima. Looks for the positive peak between the two crossing to
% negative points. It gives the value and the position. ATTENTION: The
% position is relative to the particular SO interval analysed, and not to
% the whole recording.
PosPeakSOVal = zeros(1,size(SOTemp1,2)); % Local Maxima Value Matrx
PosPeakSOLoc = zeros(1,size(SOTemp1,2)); % Local Maxima Position
% Local Minima. Same for PosPeak
NegPeakSOVal = zeros(1,size(SOTemp1,2)); % Local Minima Value Matrx
NegPeakSOLoc = zeros(1,size(SOTemp1,2));% Local Minima Position

for i = 1:size(SOTemp1,2)
    %Max
    [PosPeakSOVal(:,i), PosPeakSOLoc(:,i)] = max(FilteredEEG(SOTemp1(1,i):SOTemp1(2,i)));%,'MinPeakDistance',length((SOTemp1(1,i):SOTemp1(2,i)))/1.05); %There are cases when within the SO interval there are more than one positive or negative peak, the MinPeakDistance is set to only look for the more positive or more negative.
    %Min
    [NegPeakSOVal(:,i), NegPeakSOLoc(:,i)] = min(FilteredEEG(SOTemp1(1,i):SOTemp1(2,i)));%,'MinPeakDistance',length((SOTemp1(1,i):SOTemp1(2,i)))/1.05);
end

% Thresholds based on negpk and pk2pk
Pk2Pk       = [NegPeakSOVal; PosPeakSOVal];
Pk2PkDiff   = diff(Pk2Pk);

% Threshold independent for each Channel having as a reference Ch1
% (EEG_mPFC Left), i.e. NegPk: 36.9% and Pk2Pk: 46.8
% Distribution of the NegPks
NegPk_Edges         = (-1e-4:1e-6:0); %Values in Anuck's data are very small
[N_NegPk, E_NegPk]  = histcounts(NegPeakSOVal, NegPk_Edges);

% Distribution of the Pk2Pk
Pk2Pk_Edges         = (0:1e-6:2e-4);
[N_Pk2Pk, E_Pk2Pk]  = histcounts(Pk2PkDiff, Pk2Pk_Edges);

NegPkPerc = 36.9; % Value of events detected using mPFC_Left considering M�lle 2006 as reference.
for i = 1:length(NegPk_Edges)
    PercentCases_NegPk = sum(N_NegPk(1:i))/sum(N_NegPk)*100;
    if PercentCases_NegPk >= NegPkPerc
        ThreshNegPk = E_NegPk(i);
        PercentCases_NegPk;
        break
    end
end

Pk2PkPerc = 46.9; % Value of events detected using mPFC_Left considering M�lle 2006 as reference.
for i = 1:length(Pk2Pk_Edges)
    PercentCases_Pk2Pk = sum(N_Pk2Pk(1:i))/sum(N_Pk2Pk)*100;
    if PercentCases_Pk2Pk >= Pk2PkPerc
        ThreshPk2Pk = E_Pk2Pk(i);
        PercentCases_Pk2Pk = 100-PercentCases_Pk2Pk;
        break
    end
end

%% ValidSO fulfilling the criteria
%Creates a matrix with the beginning (1,:) and the end (2,:) of the Valid Slow Oscillations detected after Thresholds 
ValidSO = [];
j       = 1;
for i = 1:length(SOTemp1)
    if NegPeakSOVal(i) < ThreshNegPk && Pk2PkDiff(i) > ThreshPk2Pk
        ValidSO(:,j) = SOTemp1(:,i);
        j=j+1;
    end    
end

% Looking for Pos and Neg Peaks only for the Valid SO
% Local Maxima ValidSO
PosPeakValidSOVal = zeros(1,length(ValidSO)); % Local Maxima Value Matrx
PosPeakValidSOLoc = zeros(1,length(ValidSO)); % Local Maxima Position
% Local Minima ValidSO
NegPeakValidSOVal = zeros(1,length(ValidSO)); % Local Minima Value Matrx
NegPeakValidSOLoc = zeros(1,length(ValidSO)); % Local Minima Position
for i = 1:length(ValidSO)
    % Max
    [PosPeakValidSOVal(:,i), PosPeakValidSOLoc(:,i)] = max(FilteredEEG(ValidSO(1,i):ValidSO(2,i)));%,'MinPeakDistance',length((ValidSO(1,i):ValidSO(2,i)))/1.1); % To make sure that there is only 1 peak detected the MinPeakDistance is also included.
    % Min
    [NegPeakValidSOVal(:,i), NegPeakValidSOLoc(:,i)] = min(FilteredEEG(ValidSO(1,i):ValidSO(2,i)));%,'MinPeakDistance',length((ValidSO(1,i):ValidSO(2,i)))/1.1); % To make sure that there is only 1 peak detected the MinPeakDistance is also included.
end

% Valid Pk2Pk info (Amplitud)
Pk2PkValid      = [NegPeakValidSOVal; PosPeakValidSOVal];
Pk2PkValidDiff  = diff(Pk2PkValid);

%Position for the Min
ValidSONegPkLoc = zeros(size(ValidSO,2),1);
for i=1:size(ValidSO,2)
    ValidSONegPkLoc(i) = (ValidSO(1,i)+NegPeakValidSOLoc(i));
end

% Calculating the midTime - neg2pos zero crossing
j       = 1;
k       = 1;
midTime = [];
for i = 1:length(indx_up)
    if indx_up (i) > ValidSO(2,end)
        break
    end
    if indx_up(i) > ValidSO(1,j) && indx_up(i) < ValidSO(2,j)
        midTime(1,k)= indx_up(i);
        j           = j + 1;
        k           = k + 1;
    end
end

%% SO Vector
% It creates a vector with 0 when No SO are in that particular position or
% 1 when it corresponds to a SO
% It will be used to check SO relation between different Channels
SOVector = zeros(1,length(FilteredEEG));
for i = 1:length(ValidSO)
    SOVector(ValidSO(1,i):ValidSO(2,i)) = 1;
end

%% Inverting again the values for LFP signals
if strcmp(tmplabel(1:3),'CA1') == 1 || strcmp(tmplabel(1:3),'LFP') == 1
    RecordingValues = RecordingValues*-1;
end
    
%% Grand average
% It creates the matrix based on the negative peak +- 5 seconds.
SOGrandAverage      = zeros(size(ValidSO,2),2*TrialWindow*fsample);
out.trial           = {};
out.time            = {};
WindowGrandAverage  = TrialWindow*fsample; %
Cut                 = 0; %In case the Trial Window for the last detected event is longer than the recording 
CutBeg              = 0; %In case the Trial Window for the first detected starts before the recording file

for i = 1:size(ValidSO,2)
    if ((ValidSO(1,i)+NegPeakValidSOLoc(i))+WindowGrandAverage-1) > size(RecordingValues,1)
        Cut = 1;
    elseif ((ValidSO(1,i)+NegPeakValidSOLoc(i))-WindowGrandAverage-1) < 0
        CutBeg = 1;    
    else 
        SOGrandAverage(i,:)     = RecordingValues((ValidSO(1,i)+NegPeakValidSOLoc(i))-WindowGrandAverage:(ValidSO(1,i)+NegPeakValidSOLoc(i))+WindowGrandAverage-1);
        out.trial{1,i}          = SOGrandAverage(i,:)*1000;
        out.time{1,i}           = -TrialWindow:1/fsample:TrialWindow-1/fsample;
    end
end

%% SO Power
% Duration
ValidSODurationSec  = (diff(ValidSO))./fsample;

% Frequency
ValidSOFreq         = 1./ValidSODurationSec;

%% Creating eventdata
out.label           = {label};
out.fsample         = fsample;
tmptrialinfo        = zeros(size(ValidSO,2),8);
tmptrialinfo(:,1)   = (ceil(ValidSO(1,:)/1800000)*30); %To which timebin corresponds
tmptrialinfo(:,2)   = ValidSO(1,:); %startTime: begin of event (SO: first up-down zero crossing)
tmptrialinfo(:,3)   = ValidSO(2,:); %endTime: (SO: end of event second up-down zero crossing) 
tmptrialinfo(:,4)   = ValidSODurationSec; %duration: duration from start to end in seconds (SO: between the two down-to-up crossings)
tmptrialinfo(:,5)   = ValidSONegPkLoc; %minTime: time of minimum (SO: of negative half-wave/down-state*) in datapoints of original dataset
tmptrialinfo(:,6)   = NegPeakValidSOVal*1000; %minAmp: amplitude of maximum (SO: of positive half-wave/up-state*) in �V
tmptrialinfo(:,7)   = Pk2PkValidDiff*1000; %Amplitude of each event detected
tmptrialinfo(:,8)   = ValidSOFreq'; %Freq of each event detected

if Cut == 1
   tmptrialinfo = tmptrialinfo(1:size(out.time,2),:);
end

if CutBeg == 1
    out.trial   = out.trial(2:end);
    out.time    = out.time(2:end);
end

out.trialinfo       = tmptrialinfo;
out.vector          = SOVector;
out.sleeptime_sec   = ((size(find(SleepScoring == SleepStageToUse),1)) - (size(find(SleepScoring == SleepStageToUse & ScoringArtefacts == 1),1)))/fsample;
out.density_permin  = size(ValidSO,2)/(out.sleeptime_sec/60);
end