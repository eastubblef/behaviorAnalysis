
function [taskbase] = behaveLoad6workingOppStim_NewBin_LvR(fileDirectory, removeLongTrials, excludeTrials)
% 08.17.17 - trialMat has all necessary values needed - output = taskbase.trialMat to be called by other functions
% 07.25.17 - keep in mind that negative rx times are NaNs(?) & thus tossed?
% 
% 06.19.17 - updated for combining in the new rx time bins along with L v R trials: WORKS!
% 06.12.17 - updated for a new superFast bin of rxTimes < 100 ms & 100-200; 200-600; 600-1200; 1200-6000 ms
% 06.7.17  - NOTE: pXY trial 11 = Blackrock trial 12 (Blackrock = pXY + 1) for future Blackrock alignment
% 06.1.17  - updated to analyze stim trials vs. control trials in the session
% 03.06.17 - updated for rxMat as output field in taskbase structure
% 12.5.16  - updated for use with the opposite gain (starting w/ Vgatseven & eight); 
%          - NOW: a "right trial" means bar appeared to the L (and mouse moved R)
% 11.15.16 - updated for less strict incorrect trial extraction
% 11.2.16  - updated for proper trial indices for later spike-aligning
% 10.26.16 - updated for removal of trials lasting > 10s; input option 1 for yes
% 10.25.16 - updated for more reliable finding of incorrect trials
% 9.6.16   - updated for removal of any (initial) trials
% 6.01.16  - updated for removal of initial trials (ll. 147-156)
% 4.27.16  - updated for more reliable incorrect L/R trial extraction
% 4.16.16  - updated lines 397 & 398 to accomodate behavioral session 151106
% 3.28.16  - updated for initial incorrect L and R movements 
% 3.24.16  - updated for initial L and R movements that lead to reward
% 3.1.16   - use as standalone function for structure of task or to be called by recording files for alignment

% Updated BS 2.28.16 for blinking lum task & 2AFC (modified from bits of psychFitRxnTime3LumsPos3.m and Jacki's older code from July, 2015)
% Output will be only 1 structure with relevant fields "taskbase"
   
%% Initialize
posStartsTrueMaxIL = 400;
posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMinNeg = -190;
posStartsTrueMaxNeg = -200; 
posStartsTrueMinNegIR = -400;

removeLongTrials = 0;     % 1 for yes
excludeTrials = 0;        % 1 for yes
excludeWhichTrials = 186; % (currently only works from trial 1 up to this number) - work on lines 83-88 & line 99; ll. 162-172

tHold = 10;               % if velocity is desired - not in use yet
stimTypeNum = 2;             % will be 1 or 2 as of 6.1.17

fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2';
filenamestr = 'SC2';

%Define bins for rx times:
edges = [0:100:6000];
mids = conv2(edges, [1 1], 'valid')/2;
numBins = numel(edges);                  %gives 53 bins for SC2

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
binwidth = edges(2) - edges(1);          %should now be 100 ms
binwidthFastest = edges(3);              %should be 200 ms
binwidthSuperFastest = edges(2);         %100 ms

mouseName = filenamestr;
numSession = 170926;
stimEpochEnd = 500;
 
if exist('filenamestr', 'var');
    [filenamestrE, path] = uigetfile('*.csv*','select the csv file', fileDirectory);                 % get the actual # total trials 
end

if exist('filenamestr', 'var');
    [filenamestrX, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fileDirectory);          % pxy file has solenoid discharge TS with respect to trial starts
end
    
%% Load csv behavior data & implement trial structure

clear disp* trial* 

 csvFile = filenamestrE;
 data.trialTimes = dlmread(csvFile,',',2,0);

 data.attemptedTrials = data.trialTimes(:,1);                                 % gives the raw trial numbers
 data.numRewardedTrials = data.trialTimes(end,1);
 data.trialStarts = data.trialTimes(:,2);
 data.trialEnds = data.trialTimes(:,4);
 
 data.stimTime = data.trialTimes(:,7);                                        % time that the laser came on; 
 data.stimType = data.trialTimes(:,9);                                        % for stimulated trials, the stimType == 2; will register during the ITI before the stim trial starts (since it's coded to come on at trial-end + time) 
 numelStimTrialsCsv = numel(find(data.stimType == 2));                     % This should be the same number as numelStimTrialsXY (line 249)
 
 pxyFile = filenamestrX;
 data.trajectories = dlmread(pxyFile,',',2,0);                                % can use the 'taskbase' structure here to parse out solenoid TS & centered TS

 time = data.trajectories(:, 1);
 x = data.trajectories(:, 2);
 y = data.trajectories(:, 3);
 btTrials = data.trajectories(:, 5);
 luminance = data.trajectories(:, 6);                                         % for blinking bar
 stimTypeXY = data.trajectories(:,8);                                         
 
%  pFile = filenamestrF;
%  data.trialParams = dlmread(pFile,',',1+skipInitTrials,0);                  % 11th col is the column containing information about which trials were stimulated     

%% Extract relavent info: use the _pXY file

%Note, the #occurances of trial repeats not #incorrects; diff # timestamps is!
attemptTrials = data.attemptedTrials;
[Ux,xa,xc] = unique(attemptTrials,'stable');                                  % uniques
fr_x = accumarray(xc,1);                                                      % frequencies
repTrials = find(fr_x > 1);                                                   % first position of first repeating # (trials)
sumTrials = sum(fr_x);
sumReps = sumTrials - attemptTrials(end);                                     % sum of number of repeated trials (# incorrect trials)
numReps = fr_x(repTrials);                                                    % number of times repeated trials are repeated
trueNumIncorrs = sum(numReps)-numel(numReps);                                 % check if same as sumReps! Subtract 1 from ea. numReps val for true # incorrects

maxy = max(abs(y));
vel = abs(x);                                                                 % pick up any movement
%    velThold = find(vel > tHold);  NOT IN USE                                % threshold velocity, lower bound

trialStartInds = find(diff(btTrials) == -1) +1;                               % transition from 1 to 0; trial start
if excludeTrials == 1
    trialStartInds = trialStartInds(excludeWhichTrials+1:numel(trialStartInds));  % 6.1.16 for initial trial subtractions
    trialStartTimes = time(trialStartInds);
else
    trialStartTimes = time(trialStartInds);
end

% Find proper initial positions (pxy file): will be one interative value before the TS for trialStart
forPosTrialStarts = trialStartInds-1;
posTrialStarts = y(forPosTrialStarts);
%     trueTrials = find(posTrialStarts == posStartsTrueMax |...
%     posStartsTrueMaxNeg | posStartsTrueMin | posStartsTrueMinNeg);              % may need to account for premature movements before last ITI TS

if excludeTrials == 1
    trialEndInds = trialEndInds(excludeTrials+1:numel(trialEndInds));             % 6.1.16 for initial trial subtractions
else
    trialEndInds = find(diff(btTrials) == 1) +1;                                  % these are the last 0s before 1 (before next ITI)
end
trialEndTimes = time(trialEndInds);
PosTrialEnds = y(trialEndInds);
    
%% Assess ALL trials for R v L. Note, these are all before long trial subtractions

trialNum = 1:length(trialStartInds);                                          % absolute number of all trials
lumR = unique(luminance);                                             
lumL = -unique(luminance);
lumAll = vertcat(lumR, lumL);
lumTrials = luminance(trialStartInds);                                        % vector of luminance vals per trial start

%these are ALL of the initial bar/wheel starting positions
if excludeTrials == 1
    LRtrial = trialEndInds(excludeTrials+1:numel(trialEndInds));                % 6.1.16 (already inherent) for initial trial subtractions
else
    LRtrial = posTrialStarts;                                                   % absoute number of total trials
end

%% deprecated for general purposes now that recalculating usable trials
% Ltrial = posTrialStarts < 0;
% Rtrial = posTrialStarts > 0;
% numelLtrials = numel(find(LRtrial<0));
% numelRtrials = numel(find(LRtrial>0));
% 
% LtrialStarts = trialStartTimes(Ltrial == 1);
% if numel(Ltrial) > numel(trialEndTimes)
%     Ltrial = Ltrial(end-1);
% end
% LtrialEnds = trialEndTimes(Ltrial == 1);
% 
% RtrialStarts = trialStartTimes(Rtrial == 1);
% if numel(Rtrial) > numel(trialEndTimes)
%     Rtrial = Rtrial(end-1);
% end
% RtrialEnds = trialEndTimes(Rtrial == 1);
    
%% From Jacki's code to get vectors of bar position:     

position = y;                                                                   % bar position (distance for correct centering)
trialChange = diff(btTrials);                                                   % trial changes where zeros change to ones
trialchange_row = find(trialChange);                                            % give row num where trial changes

% sort trials from ITIs - find the bar position at trial start - have to go one i backward for that first value (last ITI)
barpositionCell = {};
m=1; 
% m = excludeTrials+1;
for n=1:2:(numel(trialchange_row)-1);                                           % takes only odd trials (even trials are between trials)
    barpositionCell{m} = position(trialchange_row(n):trialchange_row(n+1));     % takes CORRECTED position data for each trial and puts into indvidual cell
    time_trial{m} = time(trialchange_row(n)+1:trialchange_row(n+1));            % times for each trial - doesn't need to be corrected
    timestart{m} = time(trialchange_row(n)+1);                                  % also doesn't need correcting
    th(m) = position(trialchange_row(n));                                       % added for new "th": bar has to travel abs(th) distance to cross center 
    stimOnCell{m} = stimTypeXY(trialchange_row(n):trialchange_row(n+1));            % new - this should correct for stim coming on during previous trial end
%     stimOnCellTime{m} = time
    m=m+1;
end
    
% collect last trial (not detected by diff since no intertrial interval after last trial)
if btTrials(max(trialchange_row)+1) == 0;
    barpositionCell{m} = position(max(trialchange_row):end);                    % takes CORRECTED position data for each trial and puts into indvidual cell
    time_trial{m} = time(max(trialchange_row)+1:end);                           % times for each trial
    timestart{m} = time(max(trialchange_row)+1);
    stimOnCell{m} = stimTypeXY(max(trialchange_row):end);
end

%% Subtract initial trials for exclusion: (still need to update for stim trials)

if excludeTrials == 1 
    for e = (excludeWhichTrials)+1:numel(barpositionCell)
        newBarpositionCell{e - excludeWhichTrials} = barpositionCell{e};
        newTime_trial{e - excludeWhichTrials} = time_trial{e};
        newTimestart{e - excludeWhichTrials} = timestart{e};
    end
    barpositionCell = newBarpositionCell;
    time_trial = newTime_trial;
    timestart = newTimestart;
else
end

%% New - remove trials that lasted longer than 10s:

if removeLongTrials == 1
    for t = 1:numel(time_trial)
        trialLength(t) = time_trial{t}(end) - time_trial{t}(1);   
    end
    maxLength = max(trialLength);
    findLongTrials = trialLength > 10000;   %Now exclude them from being extracted below
    longTrials = find(findLongTrials == 1); %indices of long trials
else if removeLongTrials == 0
        maxLength = NaN;
        longTrials = NaN;
    end
end

if removeLongTrials == 1
    for z = 1:numel(barpositionCell)
       for i = 1:length(longTrials)
           if z == longTrials(i)
                barpositionCell{z} = NaN;
                time_trial{z} = NaN;
                timestart{z} = NaN;
                LRtrial(z) = NaN;
                lumTrials(z) = NaN;
           else 
           end
       end
    end
    else
end

numelLtrials = numel(find(LRtrial>0));   %new Opp task
numelRtrials = numel(find(LRtrial<0));   %new Opp task

numelRemainingTrials = numelLtrials + numelRtrials;
                     
%% determine which movement is rewarded & which trials were stimulated (6.1.17)
%   bar starts at neg. position, right WHEEL movement rewarded = 1
%   bar starts at pos. position, left WHEEL movement rewarded = -1

for n = 1:numel(barpositionCell);
    if barpositionCell{n}(1) < 0;                                               % R trials start at (-)bar positions, sorts based on first number of each cell
        reward_side(n)= -1;
    elseif barpositionCell{n}(1) > 0;                                           % L trials start at (+) bar positions, sorts based on first number of each cell
        reward_side(n)= 1;
    end
end
for p = 1:numel(stimOnCell);
    stimTrial(p) = stimOnCell{p}(1);
end
stimTrialsXY = find(stimTrial == 2);
numelStimTrialsXY = numel(find(stimTrial == 2));

%% set center to define movement as correct "choice" & find incorrects. TS are important here. 
%  correct threshold is distance from starting position to center. BUT if bar hits +/-400 (incorrect "bound"), trial is incorrect
%  BS updated to reflect 2AFC structure. n = 1xnumTrial cell array - corrected for mins (L trial) and maxes (R trial) if animal re-centers bar before solenoid discharges

% Correct R trials:
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) < 0;                                                  % R trial = L bar position & rightward movement to center
        if length(find(barpositionCell{n} >= 0) > 0) ;
            thCR(n)= min(find(barpositionCell{n} >= 0));                           % thCR = index of correct threshold crossing for R trials;
            time_thCR(n) = time_trial{n}(min(find(barpositionCell{n} >= 0))-1);    % most importantly, get the TS at thc! (have to go one backward)
        else                                                                       % if bar goes off to R on left trial
            thCR(n) = NaN;  
            time_thCR(n) = NaN;                                                    % NaN is an incorrect R trial; 0 is a L trial or one that was too long
        end                                                                        % if it's finite and > 0, its a correct R trial at that index
    end
end

% what happens to trials that I have to restart for him? - "centered" but no reward
crtNan = find(isnan(time_thCR));
crtZeros = find(time_thCR == 0);
numelCR = numel(time_thCR) - (numel(crtNan) + numel(crtZeros));
percCR = numelCR/numelRtrials * 100;

CRvalid = find(time_thCR > 1);
CRvalidRew = time_thCR(CRvalid);                                                   % centering TS for R correct

% movement extractions for correct R-movement trials
wheelCR = {};
wheelCRtime = {};
wheelRmove = {};
wheelRmoveTime = {};
wheelRfirst = [];
wheelRfirstTime = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) < 0;
        if length(find(barpositionCell{z} >= 0) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) == barpositionCell{z}(i+1) || barpositionCell{z}(i) > barpositionCell{z}(i+1)
                        wheelCR{z}(i) = NaN;
                        wheelCRtime{z}(i) = NaN;
                        wheelRmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelRmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero
                    else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                            wheelCR{z}(i) = barpositionCell{z}(i);              % wheel positions during R movements
                            wheelCRtime{z}(i) = time_trial{z}(i);               % all TS of R movements
                            wheelRmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST R movement 
                            wheelRmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of R movement
                        end
                    end
                end
            end
        end
    end
end

% initial L bar/wheel movements that lead to centering
for c = 1:length(wheelCR)
    if numel(find(wheelCR{c} >= 0) > 0)
        for i=1:length(wheelCRtime{c})
            if i == length(wheelCRtime{c})
                break
            end
            if wheelRmove{c}(i) > posStartsTrueMaxNeg && wheelRmove{c}(i) < wheelRmove{c}(i+1) && wheelRmove{c}(i) < 0 % still may need to correct for large movements that began < -200
                wheelRfirst(c) = wheelRmove{c}(i);                          % first bar position for R movements 
                wheelRfirstTime(c) = wheelRmoveTime{c}(i);                  % first TS of R movement vectors
                break                                                       % keep i from increasing once it's met this condition
            end
            
        end
    end
end

numelWheelRfirstTime = numel(find(wheelRfirstTime > 1));                    % sanity check: should be same number as numelCR var
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Incorrect R trials. For multiple min values (-400), must find the first one:  updated 10.25.16, 11.12.16
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for n=1:numel(barpositionCell);                                                    % determines position that crosses threshold
    if barpositionCell{n}(1) < 0;                                                  % R trial
        if length(find(barpositionCell{n} <= posStartsTrueMinNegIR) > 0);          % if bar goes off to R on R trial
            for i=1:length(barpositionCell{n})
                find_thIR{n}(i) = barpositionCell{n}(i);
                if length(time_trial{n}) < length(barpositionCell{n})           
                    time_trial{n}(end+1) = NaN;
                end                             
                if isnan(time_trial{n}(i+1))                                        % in case timestamp wasn't collected (160516 and test)
                    find_thIRtimes{n}(i) = time_trial{n}(i);
                    thIR(n) = numel(find_thIRtimes{n}(i));                          % index in immediately cause this is the only (or last) interesting value
                    time_thIR(n) = find_thIRtimes{n}(i);
                    n = n+1; i = 1;
                    break
                else
                   find_thIRtimes{n}(i) = time_trial{n}(i);
                end
                if find_thIR{n}(i) <= posStartsTrueMinNegIR
                    thIR(n) = find(find_thIR{n} <= posStartsTrueMinNegIR);
                    time_thIR(n) = find_thIRtimes{n}(end);
                    break
                end
                if find_thIR{n}(i) <= posStartsTrueMinNegIR && find_thIR{n}(i+1) <= posStartsTrueMinNegIR
                    find_thIR{n}(i+1) = 0;
                    find_thIRtimes{n}(i+1) = 0;
                    thIR(n) = find(find_thIR{n} <= posStartsTrueMinNegIR);
                    time_thIR(n) = find_thIRtimes{n}(i(end));
                    break
                end

            end
        else                                                                           % NaN here is a correct R trial; 0 is a L trial
            thIR(n) = NaN;      %this is the first recorded incorrect position (-400)  % if it's finite and > 0, its an incorrect R trial at that index
            time_thIR(n) = NaN; %time of first incorrect position
        end
    end
end

IRvalid = find(time_thIR > 1);     %These are the incorrect R trials         % which trials were incorrect Rs
IRvalidFail = time_thIR(IRvalid);  %And their time stamps!                   % time of bar reaching incorrect threshold

%wheel movement extractions for incorrect R trials
wheelIR = {};
wheelIRtime = {};
wheelIRmove = {};
wheelIRmoveTime = {};
wheelRfirstIR = [];
wheelRfirstTimeIR = [];
clear z;   %for troubleshooting
clear n;

for z=1:numel(barpositionCell)
    if barpositionCell{z}(1) < 0;
        if length(find(barpositionCell{z} <= posStartsTrueMinNegIR)) > 0
            for n = 1:length(find_thIR)                                      % new - correct for only 1 value found for thIR (second condition, below)
                for i=1:length(barpositionCell{z})
                    if barpositionCell{z}(i) > barpositionCell{z}(i+1) && numel(find_thIR{n} > 1)             % moving in incorrect direction
                        wheelIR{z}(i) = barpositionCell{z}(i);               % wheel positions during incorrect R trials
                        wheelIRtime{z}(i) = time_trial{z}(i);                % all TS of incorrect R trials
                        wheelIRmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST IR movement - Need to change non-position zeroes to NaNs
                        wheelIRmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of IR movement
                    end
                    if barpositionCell{z}(i) > barpositionCell{z}(i+1) && numel(find_thIR{n} == 1) % moving in incorrect direction - move this up
                        wheelIR{z}(i) = barpositionCell{z}(i);               % wheel positions during incorrect R trials
                        wheelIRtime{z}(i) = time_trial{z}(i);                % all TS of incorrect R trials
                        wheelIRmove{z}(i) = barpositionCell{z}(i);           % bit of a hack for the FIRST IR movement - Need to change non-position zeroes to NaNs
                        wheelIRmoveTime{z}(i) = time_trial{z}(i);            % & thus, hack for the FIRST TS of IR movement
                        break
                    end
                    if i+1 == length(barpositionCell{z})
                        break
                    else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                            wheelIR{z}(i) = NaN;
                            wheelIRtime{z}(i) = NaN;
                            wheelIRmove{z}(i) = NaN;                         % prevent repeating positions from being zero
                        end
                    end
                end
            end
        end
    end
end

% initial R-trial bar/wheel movements that lead to incorrect threshold
for c = 1:length(wheelIRmove)
    for i=1:length(wheelIRtime{c})
        if isfinite(wheelIRmove{c}(i)) 
            if numel(wheelIRmove{c}) == 1
                wheelRfirstIR(c) = wheelIRmove{c}(i);                             % first bar position for incorrect R movements
                wheelRfirstTimeIR(c) = wheelIRmoveTime{c}(i);                     % first TS of incorrect R movement vectors
             else if numel(wheelIRmove{c}) > 1 && wheelIRmove{c}(i) < 0
                    wheelRfirstIR(c) = wheelIRmove{c}(i);                          % first bar position for incorrect R movements
                    wheelRfirstTimeIR(c) = wheelIRmoveTime{c}(i);                  % first TS of incorrect R movement vectors
                    break                                                          % keep i from increasing once it's met this condition
                end
            end
        end
    end
end
numelWheelIRfirstTime = numel(find(wheelRfirstTimeIR > 1)); 

%% Correct L trials: 
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) > 0;                                                       % L trial
        if length(find(barpositionCell{n} <= 0) > 0) ;
            thCL(n)= min(find(abs(barpositionCell{n} <= 0)));                           % thCL = index of 1st correct threshold crossing for L trials;
            time_thCL(n) = time_trial{n}(min(find(abs(barpositionCell{n} <= 0)))-1);
        else                                                                            
            thCL(n) = NaN;              
            time_thCL(n) = NaN;                                                         % NaN is an incorrect L trial; 0 is a R trial
        end                                                                             % so if it's finite and > 0, its a correct L trial at that index
    end
end
cltNan = find(isnan(time_thCL));                                                        % sanity check
cltZeros = find(time_thCL == 0);
numelCL = numel(time_thCL) - (numel(cltNan) + numel(cltZeros));
percCL = numelCL/numelLtrials * 100;

CLvalid = find(time_thCL > 1);
CLvalidRew = time_thCL(CLvalid);                                                        % centering TSs for L correct

% L-trial movement extractions
wheelCL = {};
wheelCLtime = {};
wheelLmove = {};
wheelLmoveTime = {};
wheelLfirst = [];
wheelLfirstTime = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;
        if length(find(barpositionCell{z} <= 0) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) < barpositionCell{z}(i+1) 
                        wheelCL{z}(i) = NaN;
                        wheelCLtime{z}(i) = NaN;
                        wheelLmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelLmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero       
                    else if barpositionCell{z}(i) > barpositionCell{z}(i+1) 
                            wheelCL{z}(i) = barpositionCell{z}(i);              % wheel positions during leftward movements
                            wheelCLtime{z}(i) = time_trial{z}(i);               % all TS of leftwardward movements
                            wheelLmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST L movement 
                            wheelLmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of L movement
                        end
                    end
                end
            end
        end
    end
end

% initial L-trial movements that lead to centering 
for c = 1:length(wheelCL)
    if numel(find(wheelCL{c} <= 0) > 0)
        for i=1:length(wheelCLtime{c})
            if i == length(wheelCLtime{c})
                break
            end
            if wheelLmove{c}(i) < posStartsTrueMax && wheelLmove{c}(i) > wheelLmove{c}(i+1) && wheelLmove{c}(i) > 0 %still may need to correct for large movements that began > 200
                wheelLfirst(c) = wheelLmove{c}(i);                          % first bar position for L movements
                wheelLfirstTime(c) = wheelLmoveTime{c}(i);                  % first TS of L movement vectors
                break                                                       % keep i from increasing once it's met this condition
            end
            
        end
    end
end
numelWheelLfirstTime = numel(find(wheelLfirstTime > 1));                    % sanity check; should be same as numelCL

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Incorrect L trials:
for n=1:numel(barpositionCell);
    if barpositionCell{n}(1) > 0;                                           % L trial
        if length(find(barpositionCell{n} >= posStartsTrueMaxIL) > 0);      % if bar goes off to right on left trial
            for i=1:length(barpositionCell{n})
                find_thIL{n}(i) = barpositionCell{n}(i);
                if barpositionCell{n}(i) < 0                                % eliminate sneaky correct L trials
                    thIL(n) = NaN;
                    time_thIL(n) = NaN;
                    break
                end
                if length(time_trial{n}) < length(barpositionCell{n})       % adjust for no timestamp even though there is a position at that index
                    time_trial{n}(end+1) = NaN;
                end
                if isnan(time_trial{n}(i+1))                                        % in case timestamp wasn't collected (160516 and test)
                    find_thILtimes{n}(i) = time_trial{n}(i);
                    thIL(n) = numel(find_thILtimes{n}(i));                          % index in immediately cause this is the only (or last) interesting value
                    time_thIL(n) = find_thILtimes{n}(i);
                    n = n+1; i = 1;
                    break
                else
                    find_thILtimes{n}(i) = time_trial{n}(i);
                end
                if find_thIL{n}(i) >= posStartsTrueMaxIL
                    thIL(n) = find(find_thIL{n} >= posStartsTrueMaxIL);
                    time_thIL(n) = find_thILtimes{n}(end);
                    break
                end
                if find_thIL{n}(i) >= posStartsTrueMaxIL && find_thIL{n}(i+1) >= posStartsTrueMaxIL
                    find_thIL{n}(i+1) = 0;
                    find_thILtimes{n}(i+1) = 0;
                    thIL(n) = find(find_thIL{n} >= posStartsTrueMaxIL);
                    time_thIL(n) = find_thILtimes{n}(i(end));
                    break
                end
            end
        else                                                                    % NaN here is a correct R trial; 0 is a correct L trial
            thIL(n) = NaN;                                                      % if it's finite and > 0, its an incorrect L trial
            time_thIL(n) = NaN;
        end
    end
end
totCorr = numelCR + numelCL; % num correct trials = 250 (csv); the unaccounted-for Blacktrock trials (missing pulses) = 21; totCorr + 21 = 250

ILvalid = find(time_thIL > 1);
ILvalidFail = time_thIL(ILvalid);

%wheel movement extractions for incorrect L trials
wheelIL = {};
wheelILtime = {};
wheelILmove = {};
wheelILmoveTime = {};
wheelLfirstIL = [];
wheelLfirstTimeIL = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;                                                % So redundant... could totally refine this
        if length(find(barpositionCell{z} >= posStartsTrueMaxIL) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) > barpositionCell{z}(i+1)
                        wheelIL{z}(i) = NaN;
                        wheelILtime{z}(i) = NaN;
                        wheelILmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelILmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero
                  else if barpositionCell{z}(i) < barpositionCell{z}(i+1)        % moving in incorrect direction
                            wheelIL{z}(i) = barpositionCell{z}(i);               % wheel positions during L movements
                            wheelILtime{z}(i) = time_trial{z}(i);                % all TS of incorrect L movements
                            wheelILmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST IL movement 
                            wheelILmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of IL movement
                        end
                    end
                end
            end
        end
    end
end

% initial L bar/wheel movements that lead to incorrect threshold
clear c;
for c = 1:length(wheelILmove)
        for i=1:length(wheelILtime{c})
            if wheelILmove{c}(i) < posStartsTrueMaxIL && numel(wheelILmove{c}) < 2 % still may need to correct for large movements that began < -200
                wheelRfirstIL(c) = wheelILmove{c}(i);                          % first bar position for incorrect L movements 
                wheelRfirstTimeIL(c) = wheelILmoveTime{c}(i);                  % first TS of incorrect L movement vectors
                break                                                          % keep i from increasing once it's met this condition
            end
            if wheelILmove{c}(i) < posStartsTrueMaxIL && wheelILmove{c}(i) < wheelILmove{c}(i+1) && wheelILmove{c}(i) > 0  % still may need to correct for large movements that began < -200
                wheelLfirstIL(c) = wheelILmove{c}(i);                          % first bar position for incorrect L movements 
                wheelLfirstTimeIL(c) = wheelILmoveTime{c}(i);                  % first TS of incorrect L movement vectors
                break                                                          % keep i from increasing once it's met this condition
            end
        end
end
numelWheelILfirstTime = numel(find(wheelLfirstTimeIL > 1)); 

%% Pool all L or R movements; concatenate corrects and incorrects: 

Lcat = horzcat(wheelLfirstTime, wheelRfirstTimeIR);
allLmove = sort(Lcat);
numelLmove = numel(find(allLmove > 1));
all_Lmove = find(allLmove > 1);
all_LmoveValid = allLmove(all_Lmove);

Rcat = horzcat(wheelRfirstTime, wheelLfirstTimeIL);
allRmove = sort(Rcat);
numelRmove = numel(find(allRmove>1));
all_Rmove = find(allRmove > 1);
all_RmoveValid = allRmove(all_Rmove);

%% Tally up any trials that started & never were rewarded or incorrect (due to mouse movement and wheel centering b4 trial start)

csvTrials = data.attemptedTrials;
pxyTrials = trialStartTimes;
diffTrials = numel(pxyTrials) - numel(csvTrials);

% centering times
catCorr = [];
catCorr = horzcat(CLvalidRew, CRvalidRew);
corrects = sort(catCorr');

% bar went offscreen times:
catIncorr = [];
catIncorr = horzcat(ILvalidFail, IRvalidFail);
incorrects = sort(catIncorr');

catFailTrials1 = horzcat(ILvalid, IRvalid);
catFailTrials = sort(catFailTrials1);

numelRxTimesValid = numel(corrects) + numel(incorrects);

% post trial-removal: new calculation of trialStartInds, trialStartTimes, trialEndTimes
for t = 1:length(trialStartInds)
    if numel(trialStartInds) == numel(LRtrial) && isnan(LRtrial(t))
        newTrialStartInds(t) = NaN;
    else newTrialStartInds(t) = trialStartInds(t);
    end
    if numel(trialStartTimes) == numel(LRtrial) && isnan(LRtrial(t))
        newTrialStartTimes(t) = NaN;
    else newTrialStartTimes(t) = trialStartTimes(t);
    end
    if numel(trialEndInds) < numel(LRtrial)
        trialEndInds(end+1) = NaN;
    end
    if numel(trialEndInds) == numel(LRtrial) && isnan(LRtrial(t))
       newTrialEndInds(t) = NaN;
    else newTrialEndInds(t) = trialEndInds(t);
    end

    if numel(trialEndTimes) < numel(LRtrial)
        trialEndTimes(end+1) = NaN;
    end
    if numel(trialEndTimes) == numel(LRtrial) && isnan(LRtrial(t))
        newTrialEndTimes(t) = NaN;
    else newTrialEndTimes(t) = trialEndTimes(t);
    end
end
newTrialStartInds = newTrialStartInds';
newTrialStartTimes = newTrialStartTimes';
newTrialEndInds = newTrialEndInds';
newTrialEndTimes = newTrialEndTimes';

TSstartsNcorrects = vertcat(newTrialStartTimes, corrects);
sortTSstartsNcorrects = sort(TSstartsNcorrects);

% solenoid clicks
rewTimes = [];
newTrialEndTimesPre = isfinite(newTrialEndTimes);
newTrialEndTimes2 = newTrialEndTimes(newTrialEndTimesPre);
for t = 1:length(newTrialEndTimes2);
    for c = 1:length(corrects);                                             % solenoid click will happen bt t and t+1 for trial t; first trial is always rewarded
        if t == numel(newTrialEndTimes2)
            break
        else
            if corrects(c) > newTrialEndTimes2(t) && corrects(c) < newTrialEndTimes2(t+1)
                %             rewTimes(c-1) = trialEndTimes(t);
                rewTimes(c) = newTrialEndTimes2(t+1);
            end
        end
    end
end

rewTimesNew = find(rewTimes > 0);
rewTimes = rewTimes(rewTimesNew);
startsRewCat = horzcat(newTrialStartTimes', rewTimes);
startsRews = sort(startsRewCat);
% validStrtsRews = find(diff(startsRews)>50);                               % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
% validStrtsRews = startsRews(validStrtsRews);

numelLtrials = numel(find(LRtrial<0));
numelRtrials = numel(find(LRtrial>0));

totalNumelTrials1 = isfinite(LRtrial);
totalNumelTrials = numel(find(totalNumelTrials1 == 1));
percCorrectsTot = totCorr/totalNumelTrials * 100;                           % accounts for long trial removal

% sanity matrix: 1.trial ind, 2.trialStarts, 3.Lfirst rewarded movement, 4.Lcentered, 5.Rfirst rewarded movement, 6.Rcentered, 7.ILmoveOn, 8.barOffL, 9.IRmoveOn, 10.barOffR, 11.trial end times 12.rx times
rxTimeL = [];
rxTimeR = [];
wheelLfirstTime = wheelLfirstTime';
wheelRfirstTime = wheelRfirstTime';

trialMat = NaN(length(newTrialEndTimes), 11);                       % Fill up the matrix
if numel(newTrialEndTimes) < numel(trialStartTimes)
    trialStartTimesNew = newTrialStartTimes(1:end-1);
else trialStartTimesNew = newTrialStartTimes;
end

trialMat(:,1) = 1:numel(trialStartTimesNew);
trialMat(:,2) = trialStartTimesNew;                                 %2st col is trial start times
trialMat(:,11) = newTrialEndTimes;                                  %last col is trial end times
for i = 1:length(wheelLfirstTime)
        trialMat(i,3) = wheelLfirstTime(i);                         %3rd col is first LmoveOn
        rxTimeL(i) = wheelLfirstTime(i) - newTrialStartTimes(i);
        trialMat(i,4) = time_thCL(i);                               %4th col is time of correct T-hold crossing, Ltrials
end
for i = 1:length(wheelRfirstTime)
        trialMat(i,5) = wheelRfirstTime(i);                         %5th col is first RmoveOn
        rxTimeR(i) = wheelRfirstTime(i) - newTrialStartTimes(i);
        trialMat(i,6) = time_thCR(i);                               %6th col is time of correct T-hold crossing, Rtrials
end
for i = 1:length(wheelLfirstTimeIL)
        trialMat(i,7) = wheelLfirstTimeIL(i);                       %7th col is first ILmoveOn
        rxTimeIL(i) = wheelLfirstTimeIL(i) - newTrialStartTimes(i);
        trialMat(i,8) = time_thIL(i);                               %8th col is time of incorrect baroff, Ltrials
end
for i = 1:length(wheelRfirstTimeIR)
        trialMat(i,9) = wheelRfirstTimeIR(i);                       %9th col is first IRmoveOn
        rxTimeIR(i) = wheelRfirstTimeIR(i) - newTrialStartTimes(i);
        trialMat(i,10) = time_thIR(i);                              %10th col is time of incorrect baroff, Rtrials         
end
if numel(wheelRfirstTimeIR) == 0;
    rxTimeIR = NaN;
end

%Add an additional col at end of trialMat for was it a stim trial or not: add on line 832...
% trialMat(:,13) = stimTrial';

%These may not really be that meaningful:
rxCL = rxTimeL > -0.1;                                              %in case the first movement onset coincides with trial start for super-fast (unreal) rx times (=0)
rxCLvalid = rxTimeL(rxCL);
meanRxCLvalid = mean(rxCLvalid);
medianRxCLvalid = median(rxCLvalid);

rxCR = rxTimeR > -0.1; 
rxCRvalid = rxTimeR(rxCR);
meanRxCRvalid = mean(rxCRvalid);
medianRxCRvalid = median(rxCRvalid);

rxIL = rxTimeIL > -0.1; 
rxILvalid = rxTimeIL(rxIL);
meanRxILvalid = mean(rxILvalid);
medianRxILvalid = median(rxILvalid);

rxIR = rxTimeIR > -0.1;
rxIRvalid = rxTimeIR(rxIR);
meanRxIRvalid = mean(rxIRvalid);
medianRxIRvalid = median(rxIRvalid);
 
%% Fill a new matrix of rx times:   col 6: CL = -1; CR = 1; IL = 0; IR = 0;

rxMat = NaN(800, 7);
u = 1:totalNumelTrials;
for v = 1:length(u)                 %col1 = trial number
    rxMat(v,1) = u(v);
end
for w = 1:length(rxTimeL)
    if rxTimeL(w) < 0;              %get rid of negative rx times
        rxMat(w,2) = NaN;
    else
        rxMat(w,2) = rxTimeL(w);    %col2 = rx times for CL trials
        rxMat(w,6) = -1;
    end
end
for x = 1:length(rxTimeR)
    if rxTimeR(x) < 0;              %get rid of negative rx times
        rxMat(x,3) = NaN;
    else
        rxMat(x,3) = rxTimeR(x);    %col3 = rx times for CR trials
        rxMat(x,6) = 1;
    end
end
for y = 1:length(rxTimeIL)
    if rxTimeIL(y) < 0
        rxMat(y,4) = NaN;
    else
        rxMat(y,4) = rxTimeIL(y);    %col4 = rx times for IL trials
%         rxMat(y,6) = 0;
        rxMat(y,6) = -0.1;           %col6 ILs will now be -0.1

    end
end
for z = 1:length(rxTimeIR)
    if rxTimeIR(z) < 0
        rxMat(z,5) = NaN;
    else
        rxMat(z,5) = rxTimeIR(z);    %col5 = rx times for IR trials
%         rxMat(z,6) = 0;
        rxMat(z,6) = 0.1;            %col6 IRs will now be 0.1

    end
end
stimTrial = stimTrial';
for zStim = 1:length(stimTrial)
    if stimTrial > 3                %should never happen
        rxMat(zStim,7) = NaN;
    else 
        rxMat(zStim,7) = stimTrial(zStim);
    end                             %col 7 = whether it was a stim trial or not
end

%bin the rx times:
% hist(x,xbins)
rxCondensed = [];
for m = 1:length(rxMat(:,2))
    if isfinite(rxMat(m,2))
        rxCondensed(m,1) = rxMat(m,1);
        rxCondensed(m,2) = rxMat(m,2);
        rxCondensed(m,3) = rxMat(m,6);
        rxCondensed(m,4) = rxMat(m,7);        % for stimulated v. unstimulated trials
    else if isfinite(rxMat(m,3))
            rxCondensed(m,1) = rxMat(m,1);
            rxCondensed(m,2) = rxMat(m,3);
            rxCondensed(m,3) = rxMat(m,6);
            rxCondensed(m,4) = rxMat(m,7);
        else if isfinite(rxMat(m,4))
                rxCondensed(m,1) = rxMat(m,1);
                rxCondensed(m,2) = rxMat(m,4);
                rxCondensed(m,3) = rxMat(m,6);
                rxCondensed(m,4) = rxMat(m,7);
            else if isfinite(rxMat(m,5))
                    rxCondensed(m,1) = rxMat(m,1);
                    rxCondensed(m,2) = rxMat(m,5);
                    rxCondensed(m,3) = rxMat(m,6);
                    rxCondensed(m,4) = rxMat(m,7);
                else rxCondensed(m,1) = rxMat(m,1);
                     rxCondensed(m,2) = NaN;
                     rxCondensed(m,3) = rxMat(m,6);
                     rxCondensed(m,4) = rxMat(m,7);
                end
            end
        end
    end
end

% Sorted Rx times & drag trial info with ea:
rx = rxCondensed(:,2);
trialMat(:,12) = rx(1:length(trialMat)); %to pass Rx info for later, in-detail trial/move neural alignment

% Adjust for stimulated v. unstimulated trials:
trialMat(:,13) = stimTrial';
trialMat(:,14) = LRtrial;
trialMat(:,15) = lumTrials;

%% Pull out which wheel L/R firstTimes (correct trials) have fastest/slowest Rx times:

clear i;
rxCondensed3 = rxCondensed(:,3);
rxCondensed4 = rxCondensed(:,4);                                            %these are the stim v. unstim trials

%Add superFast rx time bin < 100 ms, L-mv:
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) < 100 && rxCondensed4(i) == 0
            wheelLfirstTimeSuperFastest(i) = wheelLfirstTime(i);                %unstim superfastest trials
        else if rxCondensed3(i) == -1 && rx(i) < 100 && rxCondensed4(i) == 2
                wheelLfirstTimeSuperFastestStim(i) = wheelLfirstTime(i);        %stim superfastest trials
            end
        end
    end
end
if exist('wheelLfirstTimeSuperFastest')
    superFastestL = wheelLfirstTimeSuperFastest > 1;
    numelSuperFastestL = numel(find(superFastestL == 1));
else
    wheelLfirstTimeSuperFastest = [];
end
if exist('wheelLfirstTimeSuperFastestStim')
    superFastestLstim = wheelLfirstTimeSuperFastestStim > 1;
    numelSuperFastestLstim = numel(find(superFastestLstim == 1));
else
    wheelLfirstTimeSuperFastestStim = [];
end

%Add superFast rx time bin < 100 ms, R-mv:
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) < 100 && rxCondensed4(j) == 0
            wheelRfirstTimeSuperFastest(j) = wheelRfirstTime(j);
        else if rxCondensed3(j) == 1 && rx(j) < 100 && rxCondensed4(j) == 2
                wheelRfirstTimeSuperFastestStim(j) = wheelRfirstTime(j);
            end
        end
    end
end
if exist('wheelRfirstTimeSuperFastest')
    superFastestR = wheelRfirstTimeSuperFastest > 1;
    numelSuperFastestR = numel(find(superFastestR == 1));
else
    wheelRfirstTimeSuperFastest = [];
end
if exist('wheelRfirstTimeSuperFastestStim')
    superFastestRstim = wheelRfirstTimeSuperFastestStim > 1;
    numelSuperFastestRstim = numel(find(superFastestRstim == 1));
else 
    wheelRfirstTimeSuperFastestStim = [];
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Updated for =>100-200 ms (fastest) time bin: L-mv
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 100 && rx(i) < 200 && rxCondensed4(i) == 0
            wheelLfirstTimeFastest(i) = wheelLfirstTime(i);
        else if rxCondensed3(i) == -1 && rx(i) >= 100 && rx(i) < 200 && rxCondensed4(i) == 2
                wheelLfirstTimeFastestStim(i) = wheelLfirstTime(i);
            end
        end
    end
end
if exist('wheelLfirstTimeFastest')
    fastestL = wheelLfirstTimeFastest > 1;
    numelFastestL = numel(find(fastestL == 1));
else
    wheelLfirstTimeFastest = [];
end
if exist('wheelLfirstTimeFastestStim')
    fastestLstim = wheelLfirstTimeFastestStim > 1;
    numelFastestLstim = numel(find(fastestLstim == 1));
else
    wheelLfirstTimeFastestStim = [];
end
   
% =>100-200 ms (fastest) time bin: R-mv 
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 100 && rx(j) < 200 && rxCondensed4(j) == 0
            wheelRfirstTimeFastest(j) = wheelRfirstTime(j);
        else if rxCondensed3(j) == 1 && rx(j) >= 100 && rx(j) < 200 && rxCondensed4(j) == 2
                wheelRfirstTimeFastestStim(j) = wheelRfirstTime(j);
            end
        end
    end
end
if exist('wheelRfirstTimeFastest')
    fastestR = wheelRfirstTimeFastest > 1;
    numelFastestR = numel(find(fastestR == 1));
else
    wheelRfirstTimeFastest = [];
end
if exist('wheelRfirstTimeFastestStim')
    fastestRstim = wheelRfirstTimeFastestStim > 1;
    numelFastestRstim = numel(find(fastestRstim == 1));
else 
    wheelRfirstTimeFastestStim = [];
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fast-rx (=>200-600ms) Left correct time starts
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 200 && rx(i) < 600 && rxCondensed4(i) == 0
            wheelLfirstTimeFast(i) = wheelLfirstTime(i);
        else if rxCondensed3(i) == -1 && rx(i) >= 200 && rx(i) < 600 && rxCondensed4(i) == 2
                wheelLfirstTimeFastStim(i) = wheelLfirstTime(i);
            end
        end
    end
end
if exist('wheelLfirstTimeFast')
    fastL = wheelLfirstTimeFast > 1;
    numelFastL = numel(find(fastL == 1));
else
    wheelLfirstTimeFast = [];
end
if exist('wheelLfirstTimeFastStim')
    fastLstim = wheelLfirstTimeFastStim > 1;
    numelFastLstim = numel(find(fastLstim == 1));
else
    wheelLfirstTimeFastStim = [];
end

%Fast (=>200-600ms) Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 200 && rx(j) < 600 && rxCondensed4(j) == 0
            wheelRfirstTimeFast(j) = wheelRfirstTime(j);
        else if rxCondensed3(j) == 1 && rx(j) >= 200 && rx(j) < 600 && rxCondensed4(j) == 2
                wheelRfirstTimeFastStim(j) = wheelRfirstTime(j);
            end
        end
    end
end
if exist('wheelRfirstTimeFast')
    fastR = wheelRfirstTimeFast > 1;
    numelFastR = numel(find(fastR == 1));
else
    wheelRfirstTimeFast = [];
end
if exist('wheelRfirstTimeFastStim')
    fastRstim = wheelRfirstTimeFastStim > 1;
    numelFastRstim = numel(find(fastRstim == 1));
else 
    wheelRfirstTimeFastStim = [];
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%mid-rx (=>600-1200)Left correct time starts
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 600 && rx(i) < 1200 && rxCondensed4(i) == 0;
            wheelLfirstTimeMids(i) = wheelLfirstTime(i);
        else if rxCondensed3(i) == -1 && rx(i) >= 600 && rx(i) < 1200 && rxCondensed4(i) == 2;
                wheelLfirstTimeMidsStim(i) = wheelLfirstTime(i);
            end
        end
    end
end

%mid-rx (=>600-1200) Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 600 && rx(j) < 1200 && rxCondensed4(j) == 0;
            wheelRfirstTimeMids(j) = wheelRfirstTime(j);
        else if rxCondensed3(j) == 1 && rx(j) >= 600 && rx(j) < 1200 && rxCondensed4(j) == 2;
                wheelRfirstTimeMidsStim(j) = wheelRfirstTime(j);
            end
        end
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%slow-rx (=>1200 - 6000)Left correct time starts
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 1200 && rx(i) < 6000 && rxCondensed4(i) == 0;
            wheelLfirstTimeSlow(i) = wheelLfirstTime(i);
        else if rxCondensed3(i) == -1 && rx(i) >= 1200 && rx(i) < 6000 && rxCondensed4(i) == 2;
                wheelLfirstTimeSlowStim(i) = wheelLfirstTime(i);
            end
        end
    end
end
if exist('wheelLfirstTimeSlow')
    slowL = wheelLfirstTimeSlow > 1;
    numelSlowL = numel(find(slowL == 1));
else wheelLfirstTimeSlow = [];
end    
if exist('wheelLfirstTimeSlowStim')
    slowLstim = wheelLfirstTimeSlowStim > 1;
    numelSlowLstim = numel(find(slowLstim == 1));
else wheelLfirstTimeSlowStim = [];
end

%slow-rx (=>1200 - 6000)Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 1200 && rx(j) < 6000 && rxCondensed4(j) == 0;
            wheelRfirstTimeSlow(j) = wheelRfirstTime(j);
        else if rxCondensed3(j) == 1 && rx(j) >= 1200 && rx(j) < 6000 && rxCondensed4(j) == 2;
                wheelRfirstTimeSlowStim(j) = wheelRfirstTime(j);
            end
        end
    end
end
if exist('wheelRfirstTimeSlow')
    slowR = wheelRfirstTimeSlow > 1;
    numelSlowR = numel(find(slowR == 1));
else wheelRfirstTimeSlow = [];
end
if exist('wheelRfirstTimeSlowStim')
    slowRstim = wheelRfirstTimeSlowStim > 1;
    numelSlowRstim = numel(find(slowRstim == 1));
else wheelRfirstTimeSlowStim = [];
end

%% Resume rx times plotting into histogram:

corrVals = rxCondensed(:,3);
[rxSort1, index] = sort(rx, 1, 'ascend');
rxSort = [];
rxSort(:,1) = index;
rxSort(:,2) = rxSort1;
rxSort(:,3) = corrVals(index);
rxSort(:,4) = rxCondensed4(index);        %add a col for the stim parameter 

% edges = [0:100:6000];                   %now defined at begining
% mids = conv2(edges, [1 1], 'valid')/2;
% numBins = numel(edges);                  %gives 53 bins for SC2

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
mu = rxSort(:,3);                          %L/R corrects; 0 = incorrects
% binwidth = edges(2) - edges(1);          %should now be 100 ms
% binwidthFastest = edges(3);              %should be 200 ms
% binwidthSuperFastest = edges(2);         %100 ms

%% Separate out the stimulated trials:
rxSort4 = rxSort(:,4);
stimMat = [];
stimMatOne = [];
for i = 1:length(rxSort(:,1))
    if rxSort(i,4) == 2
       stimMatOne(i,1) = rxSort(i,1);
    end
end
stimMatPre = stimMatOne > 1;
stimMat(:,1) = stimMatOne(stimMatPre);
stimMat(:,2) = rxSort1(stimMatPre);
stimMat(:,3) = mu(stimMatPre);
stimMat(:,4) = rxSort4(stimMatPre);

%Separate out L-mv trials within the stimulated set of trials:
stimMat1 = stimMat(:,1);
stimMat2 = stimMat(:,2);
stimMat3 = stimMat(:,3);
stimMat4 = stimMat(:,4);

stimMat_L = [];
stimMatOne_L = [];
for k = 1:length(stimMat(:,1))
    %     if stimMat(k,3) == -1                                                    %correct L-mv trials represented as -1; R-mv trials are +1
    if stimMat(k,3) < 0                                                        %ALL L-mv trials represented as < 0 now
        stimMatOne_L(k,1) = stimMat(k,1);
    end
end
stimMatPre_L = stimMatOne_L > 1;
stimMat_L(:,1) = stimMatOne_L(stimMatPre_L);
stimMat_L(:,2) = stimMat2(stimMatPre_L);
stimMat_L(:,3) = stimMat3(stimMatPre_L);
stimMat_L(:,4) = stimMat4(stimMatPre_L);

%Separate out R-mv trials within the stimulated set of trials:
stimMat_R = [];
stimMatOne_R = [];
for l = 1:length(stimMat(:,1))
    %     if stimMat(l,3) == 1                                                    %L-mv trials represented as -1; R-mv trials are +1
    if stimMat(l,3) > 0                                                    % ALL R-mv trials
        stimMatOne_R(l,1) = stimMat(l,1);
    end
end
stimMatPre_R = stimMatOne_R > 1;
stimMat_R(:,1) = stimMatOne_R(stimMatPre_R);
stimMat_R(:,2) = stimMat2(stimMatPre_R);
stimMat_R(:,3) = stimMat3(stimMatPre_R);
stimMat_R(:,4) = stimMat4(stimMatPre_R);

%% Unstimulated trials:
noStimMat = [];
noStimMatOne = [];
for j = 1:length(rxSort(:,1))
    if rxSort(j,4) == 0
       noStimMatOne(j,1) = rxSort(j,1);
    end
end
noStimMatPre = noStimMatOne > 1;
noStimMat(:,1) = noStimMatOne(noStimMatPre);
noStimMat(:,2) = rxSort1(noStimMatPre);
noStimMat(:,3) = mu(noStimMatPre);
noStimMat(:,4) = rxSort4(noStimMatPre);

%Separate out L-mv trials within the unstimulated set of trials:
noStimMat1 = noStimMat(:,1);
noStimMat2 = noStimMat(:,2);
noStimMat3 = noStimMat(:,3);
noStimMat4 = noStimMat(:,4)

noStimMat_L = [];
noStimMatOne_L = [];
for k = 1:length(noStimMat(:,1))
    %     if noStimMat(k,3) == -1                                                    %L-mv trials represented as -1; R-mv trials are +1
    if noStimMat(k,3) < 0                                                    %ALL L-mv trials represented
        noStimMatOne_L(k,1) = noStimMat(k,1);
    end
end
noStimMatPre_L = noStimMatOne_L > 1;
noStimMat_L(:,1) = noStimMatOne_L(noStimMatPre_L);
noStimMat_L(:,2) = noStimMat2(noStimMatPre_L);
noStimMat_L(:,3) = noStimMat3(noStimMatPre_L);
noStimMat_L(:,4) = noStimMat4(noStimMatPre_L);

%Separate out R-mv trials within the unstimulated set of trials:
noStimMat_R = [];
noStimMatOne_R = [];
for l = 1:length(noStimMat(:,1))
    %     if noStimMat(l,3) == 1                                                    %R-mv trials represented as -1; R-mv trials are +1
    if noStimMat(l,3) > 0                                                    %ALL R-mv trials represented
        noStimMatOne_R(l,1) = noStimMat(l,1);
    end
end
noStimMatPre_R = noStimMatOne_R > 1;
noStimMat_R(:,1) = noStimMatOne_R(noStimMatPre_R);
noStimMat_R(:,2) = noStimMat2(noStimMatPre_R);
noStimMat_R(:,3) = noStimMat3(noStimMatPre_R);
noStimMat_R(:,4) = noStimMat4(noStimMatPre_R);
 
%% Prepare the unstimulated trials for histogram-ing

noStimMat2_L = noStimMat_L(:,2);
noStimMat3_L = noStimMat_L(:,3);
noStimMat2_R = noStimMat_R(:,2);
noStimMat3_R = noStimMat_R(:,3);

%Superfastest L trials:
superFastestInds_L = find(noStimMat2_L < binwidthSuperFastest);             %superfastest trials
rxSortSuperFastest_L = noStimMat2_L(superFastestInds_L);                    %all L-mv rx times (<100ms)
% corVecSuperFastest_L = noStimMat3_L(superFastestInds_L);                   %was it a correct L mv trial
corVecSuperFastest_L1 = noStimMat3_L(superFastestInds_L);                   %All L-mv trials 
corVecSuperFastest_L = corVecSuperFastest_L1 == -1;                     
yesSuperFastest_L = find(corVecSuperFastest_L ~= 0);                        %these are superfastest L correct trials
numelCorVecSuperFastest_L = numel(corVecSuperFastest_L);                    %total num trials that were L superfastest
percSuperFastestCor_L = numel(yesSuperFastest_L)/numelCorVecSuperFastest_L; %decimal corrects for L superfastest bin
NumSuperFastest_L = percSuperFastestCor_L * 100;                            %percent corrects for L superfastest bin

%Superfastest R trials:
superFastestInds_R = find(noStimMat2_R < binwidthSuperFastest);             %superfastest R trials
rxSortSuperFastest_R = noStimMat2_R(superFastestInds_R);                    %their rx times (<100ms)
% corVecSuperFastest_R = noStimMat3_R(superFastestInds_R);
corVecSuperFastest_R1 = noStimMat3_R(superFastestInds_R);                   %All R-mv trials
corVecSuperFastest_R = corVecSuperFastest_R1 == 1;                     
yesSuperFastest_R = find(corVecSuperFastest_R ~= 0);                        %these are superfastest R correct trials
numelCorVecSuperFastest_R = numel(corVecSuperFastest_R);                    %total num trials that were superfastest
percSuperFastestCor_R = numel(yesSuperFastest_R)/numelCorVecSuperFastest_R; %decimal corrects for superfastest bin
NumSuperFastest_R = percSuperFastestCor_R * 100;                            %percent corrects for superfastest bin
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt =>100-200 ms; 
fastestInds_L = find(noStimMat2_L < binwidthFastest & noStimMat2_L >= binwidthSuperFastest);
rxSortFastest_L = noStimMat2_L(fastestInds_L);
% corVecFastest_L = noStimMat3_L(fastestInds_L);
corVecFastest_L1 = noStimMat3_L(fastestInds_L);                             %All L-mv trial
corVecFastest_L = corVecSuperFastest_L1 == -1;                     
yesFastest_L = find(corVecFastest_L ~= 0);                                  %correct L-mv trials
numelCorVecFastest_L = numel(corVecFastest_L);
percFastestCor_L = numel(yesFastest_L)/numelCorVecFastest_L;                %percent corrects for fastest bin
NumFastest_L = percFastestCor_L * 100;

fastestInds_R = find(noStimMat2_R < binwidthFastest & noStimMat2_R >= binwidthSuperFastest);
rxSortFastest_R = noStimMat2_R(fastestInds_R);
% corVecFastest_R = noStimMat3_R(fastestInds_R);
corVecFastest_R1 = noStimMat3_R(fastestInds_R);                             %All R-mv trials
corVecFastest_R = corVecSuperFastest_R1 == 1;                     
yesFastest_R = find(corVecFastest_R ~= 0);                                  %correct R-mv trials
numelCorVecFastest_R = numel(corVecFastest_R);
percFastestCor_R = numel(yesFastest_R)/numelCorVecFastest_R;                %percent corrects for fastest bin
NumFastest_R = percFastestCor_R * 100;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%In case I want just under 200 ms:

fastestInds200_L = find(noStimMat2_L <= binwidthFastest); 
rxSortFastest200_L = noStimMat2_L(fastestInds200_L);
% corVecFastest_L = noStimMat3_L(fastestInds_L);
corVecFastest200_L1 = noStimMat3_L(fastestInds200_L);                             %All L-mv trial
corVecFastest200_L = corVecFastest200_L1 == -1;                     
yesFastest200_L = find(corVecFastest200_L ~= 0);                                  %correct L-mv trials
numelCorVecFastest200_L = numel(corVecFastest200_L);
percFastestCor200_L = numel(yesFastest200_L)/numelCorVecFastest200_L;                %percent corrects for fastest bin
NumFastest200_L = percFastestCor200_L * 100;

fastestInds200_R = find(noStimMat2_R <= binwidthFastest);
rxSortFastest200_R = noStimMat2_R(fastestInds200_R);
% corVecFastest_R = noStimMat3_R(fastestInds_R);
corVecFastest200_R1 = noStimMat3_R(fastestInds200_R);                             %All R-mv trials
corVecFastest200_R = corVecFastest200_R1 == 1;                     
yesFastest200_R = find(corVecFastest200_R ~= 0);                                  %correct R-mv trials
numelCorVecFastest200_R = numel(corVecFastest200_R);
percFastestCor200_R = numel(yesFastest200_R)/numelCorVecFastest200_R;                %percent corrects for fastest bin
NumFastest200_R = percFastestCor200_R * 100;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt >=200 - 600 ms (fast bin); 
binwidthFast = edges(3):200:edges(7);                                       
fastInds_L = find(noStimMat2_L >= binwidthFast(1) & noStimMat2_L < binwidthFast(3));
rxSortFast_L = noStimMat2_L(fastInds_L);
% corVecFast_L = noStimMat3_L(fastInds_L);
corVecFast_L1 = noStimMat3_L(fastInds_L);                                   %ALL L-mv trials
corVecFast_L = corVecFast_L1 == -1;                     
yesFast_L = find(corVecFast_L ~=0);                                         %correct L-mv trials
numelCorVecFast_L = numel(corVecFast_L);
percFastCor_L = numel(yesFast_L)/numelCorVecFast_L;
NumFast_L = percFastCor_L * 100;                                            %percent corrects for fast bin

fastInds_R = find(noStimMat2_R >= binwidthFast(1) & noStimMat2_R < binwidthFast(3));
rxSortFast_R = noStimMat2_R(fastInds_R);
% corVecFast_R = noStimMat3_R(fastInds_R);
corVecFast_R1 = noStimMat3_R(fastInds_R);                                  
corVecFast_R = corVecFast_R1 == 1;                     
yesFast_R = find(corVecFast_R ~=0);
numelCorVecFast_R = numel(corVecFast_R);
percFastCor_R = numel(yesFast_R)/numelCorVecFast_R;
NumFast_R = percFastCor_R * 100;                                          
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt >=600 - 1200 ms (mid bin); 
binwidthMid = edges(7):200:edges(13);                                        
midInds_L = find(noStimMat2_L >= binwidthMid(1) & noStimMat2_L < binwidthMid(4));
rxSortMid_L = noStimMat2_L(midInds_L);
% corVecMid_L = noStimMat3_L(midInds_L);  
corVecMid_L1 = noStimMat3_L(midInds_L);                                  
corVecMid_L = corVecMid_L1 == -1;                     
yesMid_L = find(corVecMid_L ~= 0);
numelCorVecMid_L = numel(corVecMid_L);
percMidCor_L = numel(yesMid_L)/numelCorVecMid_L;
NumMid_L = percMidCor_L * 100;                              

midInds_R = find(noStimMat2_R >= binwidthMid(1) & noStimMat2_R < binwidthMid(4));
rxSortMid_R = noStimMat2_R(midInds_R);
% corVecMid_R = noStimMat3_R(midInds_R);  
corVecMid_R1 = noStimMat3_R(midInds_R);                                   
corVecMid_R = corVecMid_R1 == 1;                     
yesMid_R = find(corVecMid_R ~= 0);
numelCorVecMid_R = numel(corVecMid_R);
percMidCor_R = numel(yesMid_R)/numelCorVecMid_R;
NumMid_R = percMidCor_R * 100;                             
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

slowInds_L = find(noStimMat2_L >= binwidthMid(4));                                          
rxSortSlow_L = noStimMat2_L(slowInds_L);
% corVecSlow_L = noStimMat3_L(slowInds_L);
corVecSlow_L1 = noStimMat3_L(slowInds_L);                                  
corVecSlow_L = corVecSlow_L1 == -1;                     
yesSlow_L = find(corVecSlow_L ~= 0);
numelCorVecSlow_L = numel(corVecSlow_L);
percSlowCor_L = numel(yesSlow_L)/numelCorVecSlow_L;
NumSlow_L = percSlowCor_L * 100;                            

slowInds_R = find(noStimMat2_R >= binwidthMid(4));    
rxSortSlow_R = noStimMat2_R(slowInds_R);
% corVecSlow_R = noStimMat3_R(slowInds_R);
corVecSlow_R1 = noStimMat3_R(slowInds_R);                                   
corVecSlow_R = corVecSlow_R1 == 1;                     
yesSlow_R = find(corVecSlow_R ~= 0);
numelCorVecSlow_R = numel(corVecSlow_R);
percSlowCor_R = numel(yesSlow_R)/numelCorVecSlow_R;
NumSlow_R = percSlowCor_R * 100;     

%Trials  >=600 ms L-mv trials (mid&slow bin) for the choice matrix
binwidthMidnSlow = edges(7):200:edges(end);                                        
midnSlowInds_L = find(noStimMat2_L >= binwidthMidnSlow(1));
rxSortMidnSlow_L = noStimMat2_L(midnSlowInds_L);
% corVecMidStim_L = stimMat3_L(midIndsStim_L);  
corVecMidnSlow_L1 = noStimMat3_L(midnSlowInds_L);                        
corVecMidnSlow_L = corVecMidnSlow_L1 == -1;                     
yesMidnSlow_L = find(corVecMidnSlow_L ~= 0);
numelCorVecMidnSlow_L = numel(corVecMidnSlow_L);
percMidnSlowCor_L = numel(yesMidnSlow_L)/numelCorVecMidnSlow_L;
NumMidnSlow_L = percMidnSlowCor_L * 100;     

%Trials  >=600 ms R-mv trials (mid&slow bin) for the choice matrix
% binwidthMidnSlow = edges(7):200:edges(end);                                        
midnSlowInds_R = find(noStimMat2_R >= binwidthMidnSlow(1));
rxSortMidnSlow_R = noStimMat2_R(midnSlowInds_R);
% corVecMidStim_L = stimMat3_L(midIndsStim_L);  
corVecMidnSlow_R1 = noStimMat3_R(midnSlowInds_R);                        
corVecMidnSlow_R = corVecMidnSlow_R1 == -1;                     
yesMidnSlow_R = find(corVecMidnSlow_R ~= 0);
numelCorVecMidnSlow_R = numel(corVecMidnSlow_R);
percMidnSlowCor_R = numel(yesMidnSlow_R)/numelCorVecMidnSlow_R;
NumMidnSlow_R = percMidnSlowCor_R * 100;     

%% Do the plotting
% histRx = histc(rxSort(:,2), edges);                    % all rx times
histRx0_L = histc(rxSortSuperFastest_L, edges);              % 0-100 ms bin
histRx1_L = histc(rxSortFastest_L, edges);                   % 0-200 ms bin
histRx2_L = histc(rxSortFast_L, edges);                      % 200-600 ms bin
histRx3_L = histc(rxSortMid_L, edges);                       % 600-1200 ms bin
histRx4_L = histc(rxSortSlow_L, edges);                      % >1200 ms bin

histRx0_R = histc(rxSortSuperFastest_R, edges);              % 0-100 ms bin
histRx1_R = histc(rxSortFastest_R, edges);                   % 0-200 ms bin
histRx2_R = histc(rxSortFast_R, edges);                      % 200-600 ms bin
histRx3_R = histc(rxSortMid_R, edges);                       % 600-1200 ms bin
histRx4_R = histc(rxSortSlow_R, edges);                      % >1200 ms bin

figure; hold on; %Plot unstimulated L-mv trials
h4_L = bar(edges, histRx4_L, 'k');
h3_L = bar(edges, histRx3_L, 'g'); 
h2_L = bar(edges, histRx2_L, 'b');
h1_L = bar(edges, histRx1_L, 'r');  
h0_L = bar(edges, histRx0_L, 'y'); set(h0_L, 'BarWidth', 0.9);  % plot last for visual purposes
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [0 5000], 'Ylim', [0 50]);
    set(gca, 'Xtick', [0:500:5000]);
    num100 = num2str(round(NumSuperFastest_L)); num200 = num2str(round(NumFastest_L)); num200_600 = num2str(round(NumFast_L)); num600_1000 = num2str(round(NumMid_L)); num1000 = num2str(round(NumSlow_L));
    legend(num1000,num600_1000,num200_600,num200,num100,'Location','northeast');
    mouseInfo = strcat('all UNstim L trials', ',','% corrects', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
    hold off;

figure; hold on; %Plot unstimulated R-mv trials
h4_R = bar(edges, histRx4_R, 'k');
h3_R = bar(edges, histRx3_R, 'g'); %set(h3, 'BarWidth', 0.5);
h2_R = bar(edges, histRx2_R, 'b');
h1_R = bar(edges, histRx1_R, 'r'); 
h0_R = bar(edges, histRx0_R, 'y'); set(h0_R, 'BarWidth', 0.9);  
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [0 5000], 'Ylim', [0 50]);
    set(gca, 'Xtick', [0:500:5000]); 
    num100 = num2str(round(NumSuperFastest_R));num200 = num2str(round(NumFastest_R)); num200_600 = num2str(round(NumFast_R)); num600_1000 = num2str(round(NumMid_R)); num1000 = num2str(round(NumSlow_R));
    legend(num1000,num600_1000,num200_600,num200,num100,'Location','northeast');
    mouseInfo = strcat('all UNstim R trials', ',','% corrects', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
    hold off;

%% Prepare the stimulated trials for histogram-ing
stimMat2_L = stimMat_L(:,2);
stimMat3_L = stimMat_L(:,3);
stimMat2_R = stimMat_R(:,2);
stimMat3_R = stimMat_R(:,3);

superFastestIndsStim_L = find(stimMat2_L < binwidthSuperFastest);                   %superfastest trials
rxSortSuperFastestStim_L = stimMat2_L(superFastestIndsStim_L);                      %their rx times (<100ms)
% corVecSuperFastestStim_L = stimMat3_L(superFastestIndsStim_L);                      
corVecSuperFastestStim_L1 = stimMat3_L(superFastestIndsStim_L);                       
corVecSuperFastestStim_L = corVecSuperFastestStim_L1 == -1;                     
yesSuperFastestStim_L = find(corVecSuperFastestStim_L ~= 0);                        %these are superfastest correct trials
numelCorVecSuperFastestStim_L = numel(corVecSuperFastestStim_L);                    %total num trials that were superfastest
percSuperFastestCorStim_L = numel(yesSuperFastestStim_L)/numelCorVecSuperFastestStim_L; %decimal corrects for superfastest bin
NumSuperFastestStim_L = percSuperFastestCorStim_L * 100;                            %percent corrects for superfastest bin

superFastestIndsStim_R = find(stimMat2_R < binwidthSuperFastest);                   %superfastest trials
rxSortSuperFastestStim_R = stimMat2_R(superFastestIndsStim_R);                      %their rx times (<100ms)
% corVecSuperFastestStim_R = stimMat3_R(superFastestIndsStim_R);                    %was it an R mv trial
corVecSuperFastestStim_R1 = stimMat3_R(superFastestIndsStim_R);                     %was it a correct R mv trial
corVecSuperFastestStim_R = corVecSuperFastestStim_R1 == 1;                     
yesSuperFastestStim_R = find(corVecSuperFastestStim_R ~= 0);                        %these are superfastest correct trials
numelCorVecSuperFastestStim_R = numel(corVecSuperFastestStim_R);                    %total num trials that were superfastest
percSuperFastestCorStim_R = numel(yesSuperFastestStim_R)/numelCorVecSuperFastestStim_R; %decimal corrects for superfastest bin
NumSuperFastestStim_R = percSuperFastestCorStim_R * 100;                            %percent corrects for superfastest bin
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt =>100-200 ms; 
fastestIndsStim_L = find(stimMat2_L < binwidthFastest & stimMat2_L >= binwidthSuperFastest);
rxSortFastestStim_L = stimMat2_L(fastestIndsStim_L);
% corVecFastestStim_L = stimMat3_L(fastestIndsStim_L);
corVecFastestStim_L1 = stimMat3_L(fastestIndsStim_L);                        
corVecFastestStim_L = corVecFastestStim_L1 == -1;                     
yesFastestStim_L = find(corVecFastestStim_L ~= 0); 
numelCorVecFastestStim_L = numel(corVecFastestStim_L);
percFastestCorStim_L = numel(yesFastestStim_L)/numelCorVecFastestStim_L;   
NumFastestStim_L = percFastestCorStim_L * 100;

fastestIndsStim_R = find(stimMat2_R < binwidthFastest & stimMat2_R >= binwidthSuperFastest);
rxSortFastestStim_R = stimMat2_R(fastestIndsStim_R);
% corVecFastestStim_R = stimMat3_R(fastestIndsStim_R);
corVecFastestStim_R1 = stimMat3_R(superFastestIndsStim_R);                         
corVecFastestStim_R = corVecFastestStim_R1 == 1;                     
yesFastestStim_R = find(corVecFastestStim_R ~= 0); 
numelCorVecFastestStim_R = numel(corVecFastestStim_R);
percFastestCorStim_R = numel(yesFastestStim_R)/numelCorVecFastestStim_R;  
NumFastestStim_R = percFastestCorStim_R * 100;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials under 200 ms; 
fastestIndsStim200_L = find(stimMat2_L < binwidthFastest);
rxSortFastestStim200_L = stimMat2_L(fastestIndsStim200_L);
% corVecFastestStim_L = stimMat3_L(fastestIndsStim_L);
corVecFastestStim200_L1 = stimMat3_L(fastestIndsStim200_L);                        
corVecFastestStim200_L = corVecFastestStim200_L1 == -1;                     
yesFastestStim200_L = find(corVecFastestStim200_L ~= 0); 
numelCorVecFastestStim200_L = numel(corVecFastestStim200_L);
percFastestCorStim200_L = numel(yesFastestStim200_L)/numelCorVecFastestStim200_L;   
NumFastestStim200_L = percFastestCorStim200_L * 100;

fastestIndsStim200_R = find(stimMat2_R < binwidthFastest);
rxSortFastestStim200_R = stimMat2_R(fastestIndsStim200_R);
% corVecFastestStim_R = stimMat3_R(fastestIndsStim_R);
corVecFastestStim200_R1 = stimMat3_R(fastestIndsStim200_R);                         
corVecFastestStim200_R = corVecFastestStim200_R1 == 1;                     
yesFastestStim200_R = find(corVecFastestStim200_R ~= 0); 
numelCorVecFastestStim200_R = numel(corVecFastestStim200_R);
percFastestCorStim200_R = numel(yesFastestStim200_R)/numelCorVecFastestStim200_R;  
NumFastestStim200_R = percFastestCorStim200_R * 100;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt >=200 - 600 ms (fast bin); 
binwidthFast = edges(3):200:edges(7);                                    
fastIndsStim_L = find(stimMat2_L >= binwidthFast(1) & stimMat2_L < binwidthFast(3));
rxSortFastStim_L = stimMat2_L(fastIndsStim_L);
% corVecFastStim_L = stimMat3_L(fastIndsStim_L);
corVecFastStim_L1 = stimMat3_L(fastIndsStim_L);                         
corVecFastStim_L = corVecFastStim_L1 == -1;                     
yesFastStim_L = find(corVecFastStim_L ~=0);
numelCorVecFastStim_L = numel(corVecFastStim_L);
percFastCorStim_L = numel(yesFastStim_L)/numelCorVecFastStim_L;
NumFastStim_L = percFastCorStim_L * 100;                             

fastIndsStim_R = find(stimMat2_R >= binwidthFast(1) & stimMat2_R < binwidthFast(3));
rxSortFastStim_R = stimMat2_R(fastIndsStim_R);
% corVecFastStim_R = stimMat3_R(fastIndsStim_R);
corVecFastStim_R1 = stimMat3_R(fastIndsStim_R);                         
corVecFastStim_R = corVecFastStim_R1 == 1;                     
yesFastStim_R = find(corVecFastStim_R ~=0);
numelCorVecFastStim_R = numel(corVecFastStim_R);
percFastCorStim_R = numel(yesFastStim_R)/numelCorVecFastStim_R;
NumFastStim_R = percFastCorStim_R * 100;                                   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Trials bt >=600 - 1200 ms (mid bin); 
binwidthMid = edges(7):200:edges(13);                                        
midIndsStim_L = find(stimMat2_L >= binwidthMid(1) & stimMat2_L < binwidthMid(4));
rxSortMidStim_L = stimMat2_L(midIndsStim_L);
% corVecMidStim_L = stimMat3_L(midIndsStim_L);  
corVecMidStim_L1 = stimMat3_L(midIndsStim_L);                        
corVecMidStim_L = corVecMidStim_L1 == -1;                     
yesMidStim_L = find(corVecMidStim_L ~= 0);
numelCorVecMidStim_L = numel(corVecMidStim_L);
percMidCorStim_L = numel(yesMidStim_L)/numelCorVecMidStim_L;
NumMidStim_L = percMidCorStim_L * 100;                              

midIndsStim_R = find(stimMat2_R >= binwidthMid(1) & stimMat2_R < binwidthMid(4));
rxSortMidStim_R = stimMat2_R(midIndsStim_R);
% corVecMidStim_R = stimMat3_R(midIndsStim_R);  
corVecMidStim_R1 = stimMat3_R(midIndsStim_R);                        
corVecMidStim_R = corVecMidStim_R1 == 1;                     
yesMidStim_R = find(corVecMidStim_R ~= 0);
numelCorVecMidStim_R = numel(corVecMidStim_R);
percMidCorStim_R = numel(yesMidStim_R)/numelCorVecMidStim_R;
NumMidStim_R = percMidCorStim_R * 100;                              
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

slowIndsStim_L = find(stimMat2_L >= binwidthMid(4));
rxSortSlowStim_L = stimMat2_L(slowIndsStim_L);
% corVecSlowStim_L = stimMat3_L(slowIndsStim_L);
corVecSlowStim_L1 = stimMat3_L(slowIndsStim_L);                        
corVecSlowStim_L = corVecSlowStim_L1 == -1;                     
yesSlowStim_L = find(corVecSlowStim_L ~= 0);
numelCorVecSlowStim_L = numel(corVecSlowStim_L);
percSlowCorStim_L = numel(yesSlowStim_L)/numelCorVecSlowStim_L;
NumSlowStim_L = percSlowCorStim_L * 100;                            

slowIndsStim_R = find(stimMat2_R >= binwidthMid(4));
rxSortSlowStim_R = stimMat2_R(slowIndsStim_R);
% corVecSlowStim_R = stimMat3_R(slowIndsStim_R);
corVecSlowStim_R1 = stimMat3_R(slowIndsStim_R);                        
corVecSlowStim_R = corVecSlowStim_R1 == 1;                     
yesSlowStim_R = find(corVecSlowStim_R ~= 0);
numelCorVecSlowStim_R = numel(corVecSlowStim_R);
percSlowCorStim_R = numel(yesSlowStim_R)/numelCorVecSlowStim_R;
NumSlowStim_R = percSlowCorStim_R * 100;     

%Trials  >=600 ms L-mv trials (mid&slow bin) for the choice matrix
binwidthMidnSlow = edges(7):200:edges(end);                                        
midnSlowIndsStim_L = find(stimMat2_L >= binwidthMidnSlow(1));
rxSortMidnSlowStim_L = stimMat2_L(midnSlowIndsStim_L);
% corVecMidStim_L = stimMat3_L(midIndsStim_L);  
corVecMidnSlowStim_L1 = stimMat3_L(midnSlowIndsStim_L);                        
corVecMidnSlowStim_L = corVecMidnSlowStim_L1 == -1;                     
yesMidnSlowStim_L = find(corVecMidnSlowStim_L ~= 0);
numelCorVecMidnSlowStim_L = numel(corVecMidnSlowStim_L);
percMidnSlowCorStim_L = numel(yesMidnSlowStim_L)/numelCorVecMidnSlowStim_L;
NumMidnSlowStim_L = percMidnSlowCorStim_L * 100;     

%Trials  >=600 ms R-mv trials (mid&slow bin) for the choice matrix
midnSlowIndsStim_R = find(stimMat2_R >= binwidthMidnSlow(1));
rxSortMidnSlowStim_R = stimMat2_R(midnSlowIndsStim_R);
% corVecMidStim_L = stimMat3_L(midIndsStim_L);  
corVecMidnSlowStim_R1 = stimMat3_R(midnSlowIndsStim_R);                        
corVecMidnSlowStim_R = corVecMidnSlowStim_R1 == -1;                     
yesMidnSlowStim_R = find(corVecMidnSlowStim_R ~= 0);
numelCorVecMidnSlowStim_R = numel(corVecMidnSlowStim_R);
percMidnSlowCorStim_R = numel(yesMidnSlowStim_R)/numelCorVecMidnSlowStim_R;
NumMidnSlowStim_R = percMidnSlowCorStim_R * 100;    

%Do the plotting for stim trials:
histRx0stim_L = histc(rxSortSuperFastestStim_L, edges);              % 0-100 ms bin
histRx1stim_L = histc(rxSortFastestStim_L, edges);                   % 100-200 ms bin
histRx2stim_L = histc(rxSortFastStim_L, edges);                      % 200-600 ms bin
histRx3stim_L = histc(rxSortMidStim_L, edges);                       % 600-1200 ms bin
histRx4stim_L = histc(rxSortSlowStim_L, edges);                      % >1200 ms bin

histRx0stim_R = histc(rxSortSuperFastestStim_R, edges);              % 0-100 ms bin
histRx1stim_R = histc(rxSortFastestStim_R, edges);                   % 100-200 ms bin
histRx2stim_R = histc(rxSortFastStim_R, edges);                      % 200-600 ms bin
histRx3stim_R = histc(rxSortMidStim_R, edges);                       % 600-1200 ms bin
histRx4stim_R = histc(rxSortSlowStim_R, edges);                      % >1200 ms bin

figure; hold on;
h4stim_L = bar(edges, histRx4stim_L, 'k');
h3stim_L = bar(edges, histRx3stim_L, 'g'); 
h2stim_L = bar(edges, histRx2stim_L, 'b');
h1stim_L = bar(edges, histRx1stim_L, 'r');  
h0stim_L = bar(edges, histRx0stim_L, 'y'); set(h0stim_L, 'BarWidth', 0.9);  
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [0 5000], 'Ylim', [0 50]);
    set(gca, 'Xtick', [0:500:5000]); 
    num100 = num2str(round(NumSuperFastestStim_L));num200 = num2str(round(NumFastestStim_L)); num200_600 = num2str(round(NumFastStim_L)); num600_1000 = num2str(round(NumMidStim_L)); num1000 = num2str(round(NumSlowStim_L));
    legend(num1000,num600_1000,num200_600,num200,num100,'Location','northeast');
    mouseInfo = strcat('all L STIM trials', ',','% corrects', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
    hold off;
 
figure; hold on;
h4stim_R = bar(edges, histRx4stim_R, 'k');
h3stim_R = bar(edges, histRx3stim_R, 'g'); 
h2stim_R = bar(edges, histRx2stim_R, 'b');
h1stim_R = bar(edges, histRx1stim_R, 'r'); set(h1stim_R, 'BarWidth', 0.8);  
h0stim_R = bar(edges, histRx0stim_R, 'y'); set(h0stim_R, 'BarWidth',0.8);  
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [0 5000], 'Ylim', [0 50]);
    set(gca, 'Xtick', [0:500:5000]); 
    num100 = num2str(round(NumSuperFastestStim_R));num200 = num2str(round(NumFastestStim_R)); num200_600 = num2str(round(NumFastStim_R)); num600_1200 = num2str(round(NumMidStim_R)); num1200 = num2str(round(NumSlowStim_R));
    legend(num1200,num600_1200,num200_600,num200,num100,'Location','northeast');
    mouseInfo = strcat('all R STIM trials', ',','% corrects', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);

    
%% Plot the steadiness of correct trials over the session:

% figure(); plot(cumsum(diff(unique([CRvalid CLvalid])))); hold on; plot(diff(cumsum([1:537])));
% percCorrJosh = numel(unique([CRvalid CLvalid]))./537

%% Run a Wilcoxon rank sum (sign rank?) test for the stim R v. unstim Rx distributions
% [pRank_R,hRank_R,statsRank_R] = signrank(stimMat2_R, noStimMat2_R);      %nope, x & y must be same size
% [pRank_L,hRank_L,statsRank_L] = signrank(stimMat2_L, noStimMat2_L);

% [pDist,hDist,statsDist] = ranksum(histRxUnstim,histRxStim);
[pDist_L,hDist_L,statsDist_L] = ranksum(stimMat2_L, noStimMat2_L);         %these include correct & incorrect trials

figure; hold on; 
[NstimL,EstimL] = histcounts(stimMat2_L)
EstimL = EstimL(1:end-1);
[NnoStimL,EnoStimL] = histcounts(noStimMat2_L)
EnoStimL = EnoStimL(1:end-1);
plot(EstimL,NstimL, 'bl-');
plot(EnoStimL,NnoStimL, 'k-');
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [-10 3000], 'Ylim', [0 100]);
    Info = strcat('all L trials', ',','rx times', ',', mouseName, ',', num2str(numSession));
    title(Info);
figure; hold on;
[NstimR,EstimR] = histcounts(stimMat2_R)
EstimR = EstimR(1:end-1);
[NnoStimR,EnoStimR] = histcounts(noStimMat2_R)
EnoStimR = EnoStimR(1:end-1);
plot(EstimR,NstimR, 'bl-');
plot(EnoStimR,NnoStimR, 'k-');
    xlabel('Rx times (ms)'); ylabel('# Trials');
    set(gca, 'Xlim', [-10 3000], 'Ylim', [0 100]);
    Info = strcat('all R trials', ',','rx times', ',', mouseName, ',', num2str(numSession));
    title(Info);

[pDist_R,hDist_R,statsDist_R] = ranksum(stimMat2_R, noStimMat2_R);         %these include correct & incorrect trials

%% Run a test for choice on the < 200 ms rx times stim v unstim; same for 200-600 ms:
% [pDist_200L,hDist_200L,statsDist_200L] = ranksum(corVecFastest200_L, corVecFastestStim200_L);
% [pDist_200R,hDist_200R,statsDist_200R] = ranksum(corVecFastest200_R, corVecFastestStim200_R);
% 
% [pDist_200_600L,hDist_200_600L,statsDist_200_600L] = ranksum(corVecFast_L, corVecFastStim_L);
% [pDist_200_600R,hDist_200_600R,statsDist_200_600R] = ranksum(corVecFast_R, corVecFastStim_R);

%% Create a choice matrix: col1 = should-mv L trials; col2 = should-mv R trials; row 1 = unstim %corr; row 2 = stim % corr
%for the 200-600 ms rx times:
choiceMatSlow(1,1) = numelCorVecMidnSlow_L;
choiceMatSlow(1,2) = numelCorVecMidnSlow_R;
choiceMatSlow(2,2) = percMidnSlowCor_L;
choiceMatSlow(2,2) = percMidnSlowCor_R;
choiceMatSlow(3,1) = numelCorVecMidnSlowStim_L;
choiceMatSlow(3,2) = numelCorVecMidnSlowStim_R;
choiceMatSlow(4,1) = percMidnSlowCorStim_L;
choiceMatSlow(4,2) = percMidnSlowCorStim_R;

choiceMat2to600s(1,1) = numelCorVecFast_L;          %numel delivered; unstim L 
choiceMat2to600s(1,2) = numelCorVecFast_R;          %numel delivered; unstim R 
choiceMat2to600s(2,1) = percFastCor_L;              %unstim L
choiceMat2to600s(2,2) = percFastCor_R;              %unstim R
choiceMat2to600s(3,1) = numelCorVecFastStim_L;
choiceMat2to600s(3,2) = numelCorVecFastStim_R;
choiceMat2to600s(4,1) = percFastCorStim_L;          %stim L
choiceMat2to600s(4,2) = percFastCorStim_R;          %stim R

%for the <= 200 ms rx times:
choiceMat200s(1,1) = numelCorVecFastest200_L;       %numel delivered; unstim L 
choiceMat200s(1,2) = numelCorVecFastest200_R;       %numel delivered; unstim R
choiceMat200s(2,1) = percFastestCor200_L;           %unstim L
choiceMat200s(2,2) = percFastestCor200_R;           %unstim R
choiceMat200s(3,1) = numelCorVecFastestStim200_L;
choiceMat200s(3,2) = numelCorVecFastestStim200_R;
choiceMat200s(4,1) = percFastestCorStim200_L;       %stim L
choiceMat200s(4,2) = percFastestCorStim200_R;       %stim R

%% Create an rx-time fig for Gidon of actual L v R movements (independent of what kind of trial it was)

stimWentR = stimMat3 == 1 | stimMat3 == -0.1;
stimWentRrx = stimMat2(stimWentR);    %this vector

stimWentL = stimMat3 == -1 | stimMat3 == 0.1;
stimWentLrx = stimMat2(stimWentL);    %this vector
% NaNStim = stimMat2 == NaN;
% numelNaNStim = numel(NaNStim == 1);
numelNotNanStim = numel(stimWentRrx) + numel(stimWentLrx); %get rid of the NaN trials
% fracStimWentRrx = numel(stimWentRrx)/numelNotNanStim; %fraction of trials out of all stimulated that they went R
% fracStimWentLrx = numel(stimWentLrx)/numelNotNanStim; %fraction of trials out of all stimulated that they went L

noStimWentR = noStimMat3 == 1 | noStimMat3 == -0.1;
noStimWentRrx = noStimMat2(noStimWentR);  %this vector

noStimWentL = noStimMat3 == -1 | noStimMat3 == 0.1; 
noStimWentLrx = noStimMat2(noStimWentL);  %this vector
numelNotNanNoStim = numel(noStimWentRrx) + numel(noStimWentLrx);

% fracNoStimWentRrx = numel(noStimWentRrx)/numelNotNanNoStim; %fraction of trials out of unstimulated that they went R
% fracNoStimWentLrx = numel(noStimWentLrx)/numelNotNanNoStim; %fraction of trials out of all stimulated that they went L
numelNotNanWentR = numel(stimWentRrx) + numel(noStimWentRrx);
numelNotNanWentL = numel(stimWentLrx) + numel(noStimWentLrx);
fracWentRStim = numel(stimWentRrx)/numelNotNanWentR;
fracWentRnoStim = numel(noStimWentRrx)/numelNotNanWentR;
fracWentLStim = numel(stimWentLrx)/numelNotNanWentL;
fracWentLnoStim = numel(noStimWentLrx)/numelNotNanWentL;

[pDistChoice_L] = ranksum(noStimWentLrx, stimWentLrx);         %these include correct & incorrect trials
[pDistChoice_R] = ranksum(noStimWentRrx, stimWentRrx);         %these include correct & incorrect trials

figure; hold on;
plot(1, noStimWentLrx,'k.', 'markersize', 10);
plot(2, stimWentLrx, 'bl.', 'markersize', 10);
% xticklabel(['no stim', 'stim']); 
ylabel('Rx times (ms) went L');
    set(gca, 'Xlim', [0.5 2.5], 'Ylim', [0 3000]);
    set(gca, 'Xtick', [1 2]); 
    mouseInfo = strcat('trials mouse went L stim v unstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);

figure; hold on;
plot(1, noStimWentRrx,'k.', 'markersize', 10);
plot(2, stimWentRrx, 'bl.', 'markersize', 10);
% xlabel(['no stim', 'stim']); 
ylabel('Rx times (ms) went R');
    set(gca, 'Xlim', [0.5 2.5], 'Ylim', [0 3000]);
    set(gca, 'Xtick', [1 2]); 
    mouseInfo = strcat('trials mouse went R stim v unstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);

%Tease out only the trials that occurred w/in the stim epoch:
stimEndRxLnoStim1 = find(noStimWentLrx < stimEpochEnd);
stimEndRxLnoStim = noStimWentLrx(stimEndRxLnoStim1)

stimEndRxLstim1 = find(stimWentLrx < stimEpochEnd);
stimEndRxLstim = stimWentLrx(stimEndRxLstim1);

stimEndRxRnoStim1 = find(noStimWentRrx < stimEpochEnd);
stimEndRxRnoStim = noStimWentRrx(stimEndRxRnoStim1);

stimEndRxRstim1 = find(stimWentRrx < stimEpochEnd);
stimEndRxRstim = stimWentRrx(stimEndRxRstim1);
    
[pDistChoiceStimEnd_L] = ranksum(stimEndRxLnoStim, stimEndRxLstim);         %these include correct & incorrect trials
[pDistChoiceStimEnd_R] = ranksum(stimEndRxRnoStim, stimEndRxRstim);         %these include correct & incorrect trials
   
figure; hold on;
plot(1, stimEndRxLnoStim,'k.', 'markersize', 10);
plot(2, stimEndRxLstim, 'bl.', 'markersize', 10);
% xticklabel(['no stim', 'stim']); 
ylabel('Rx times (ms) went L');
    set(gca, 'Xlim', [0.5 2.5], 'Ylim', [0 1500]);
    set(gca, 'Xtick', [1 2]); 
    mouseInfo = strcat('trials mouse went L < stimEpoch; stim v unstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);

figure; hold on;
plot(1, stimEndRxRnoStim,'k.', 'markersize', 10);
plot(2, stimEndRxRstim, 'bl.', 'markersize', 10);
% xlabel(['no stim', 'stim']); 
ylabel('Rx times (ms) went R');
    set(gca, 'Xlim', [0.5 2.5], 'Ylim', [0 1500]);
    set(gca, 'Xtick', [1 2]); 
    mouseInfo = strcat('trials mouse went R < stimEpoch; stim v unstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);



%% Create taskbase structure (row vectors)

% taskbase.rewTimes = rewTimes(2:end);     % correct for zeros at first index - fixed at line 694

taskbase.rewTimes = rewTimes;
taskbase.csvTrialInds = csvTrials';
% taskbase.LtrialStarts = LtrialStarts'; % no longer accurate
% taskbase.RtrialStarts = RtrialStarts'; % no longer accurate
taskbase.leftCenteredTimes = CLvalidRew;
taskbase.leftIndsCenteredTimes = time_thCL;
taskbase.leftIncorrectTimes = ILvalidFail;
taskbase.rightCenteredTimes = CRvalidRew;
taskbase.rightIndsCenteredTimes = time_thCR;
taskbase.rightIncorrectTimes = IRvalidFail;
taskbase.trialStartInds = newTrialStartInds;
taskbase.trialStartTimes = newTrialStartTimes';
taskbase.diffTrials = diffTrials;                                           % if below number doesn't = num Blackrock pulses, add this num to it
taskbase.events = sortTSstartsNcorrects';                                   % numel here should = num pulses into blackrock for .nev file
taskbase.startsRews = startsRews;
taskbase.LRtrial = LRtrial';
taskbase.wheelLfirstTime = wheelLfirstTime';
taskbase.wheelLfirst = wheelLfirst;
taskbase.wheelILfirst = wheelLfirstIL;
taskbase.wheelILfirstTime = wheelLfirstTimeIL;
taskbase.wheelRfirst = wheelRfirst;
taskbase.wheelRfirstTime = wheelRfirstTime';
taskbase.wheelIRfirst = wheelRfirstIR;
taskbase.wheelIRfirstTime = wheelRfirstTimeIR;
taskbase.wheelAllLeft = all_LmoveValid;
taskbase.wheelAllRight = all_RmoveValid;
taskbase.trialEndTimes = newTrialEndTimes2;
taskbase.trialMat = trialMat;

taskbase.rxMat = rxMat;
taskbase.rxCondensed = rxCondensed;
taskbase.stimMat = stimMat;
taskbase.noStimMat = noStimMat;

taskbase.wheelLfirstTimeFastest = wheelLfirstTimeFastest;
taskbase.wheelRfirstTimeFastest = wheelRfirstTimeFastest;
% taskbase.wheelLfirstTimeMids = wheelLfirstTimeMids;
% taskbase.wheelRfirstTimeMids = wheelRfirstTimeMids;
taskbase.wheelLfirstTimeSlow = wheelLfirstTimeSlow;
taskbase.wheelRfirstTimeSlow = wheelRfirstTimeSlow;

taskbase.barpositionCell = barpositionCell;
taskbase.time_trial = time_trial;
taskbase.wheelLmove = wheelLmove;
taskbase.wheelRmove = wheelRmove;
taskbase.wheelILmove = wheelILmove;
taskbase.wheelIRmove = wheelIRmove;


%% COMPLETED ALL ANALYSIS. SAVE taskase STRUCTURE

disp(' ');    
disp('%-------------------------------------------------------------------');
disp(['Completed file: ' filenamestrX]);
  
filenamestrNew = fileDirectory;
targetName = filenamestrX(1:end-11);     
% cd(filenamestrNew);
     
 save([targetName '_tb.mat'],'taskbase');                              
 disp(['saved as ' targetName '_tb.mat']);
                                                                
%     save(FR_file{n_i}, 'FR', 'FR_bine');
%     save([targetName '_bh.mat'],'ContData');      
     
disp('%-------------------------------------------------------------------');
disp(' ');
    
