
function [taskbase] = behaveLoad6working(fileDirectory, removeLongTrials, excludeTrials)

% 03.06.17 - updated for rxMat as output field in taskbase structure
% 12.20.16 - updated for more meaningful trialMat and removal of reaction times that are too fast (near the end)
% 11.15.16 - updated for more reliable incorrect trial extraction
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
posStartsTrueMaxIR = 400;
posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMinNeg = -190;
posStartsTrueMaxNeg = -200; 
posStartsTrueMinNegIL = -400;

% posStartsTrueMaxIR = 350;
% posStartsTrueMax = 180;
% posStartsTrueMin = 170;
% posStartsTrueMinNeg = -170;
% posStartsTrueMaxNeg = -180; 
% posStartsTrueMinNegIL = -350;

removeLongTrials = 0;     % 1 for yes
excludeTrials = 0;        % 1 for yes
excludeWhichTrials = 186; % (currently only works from trial 1 up to this number) - work on lines 83-88 & line 99; ll. 162-172

tHold = 10;     % if velocity is desired - not in use yet

% fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/Vgatfour/forAnalysis';
% filenamestr = 'Vgatfour';
% fileDirectory = '/Users/stubblefielde/Desktop/behaviorOctNov15/Vgatthree';
% filenamestr = 'Vgatthree';
fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/Vgatfive/forAnalysis;';
filenamestr = 'Vgatfive';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170111';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings16_17/170112/1stPass/behaveSegs1';
% filenamestr  = 'Vgattwotag';
% % fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave/untag/newBehaveUnits/untagged';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/behave';
% filenamestr  = 'Vgattwotag';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks/newBehaveUnits/behaveChunks';
% filenamestr  = 'Vgatthreetag';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/newBehaveUnits/behaveChunks;'
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behave/newUnits';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104';
mouseName = filenamestr;
numSession = 151106;

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

 pxyFile = filenamestrX;
 data.trajectories = dlmread(pxyFile,',',2,0);                                % can use the 'taskbase' structure here to parse out solenoid TS & centered TS

 time = data.trajectories(:, 1);
 x = data.trajectories(:, 2);
 y = data.trajectories(:, 3);
 btTrials = data.trajectories(:, 5);
 luminance = data.trajectories(:, 7);                                         % for blinking bar

%  pFile = filenamestrF;
%  data.trialParams = dlmread(pFile,',',1+skipInitTrials,0);                  % 11th col is the column containing information about which trials were stimulated     

%% Extract relavent info: use the _pXY file

%Note, the #occurances of trial repeats not #incorrects; diff # timestamps is!
% attemptTrials = data.attemptedTrials;
% [Ux,xa,xc] = unique(attemptTrials,'stable');                                  % uniques
% fr_x = accumarray(xc,1);                                                      % frequencies
% repTrials = find(fr_x > 1);                                                   % first position of first repeating # (trials)
% sumTrials = sum(fr_x);
% sumReps = sumTrials - attemptTrials(end);                                     % sum of number of repeated trials (# incorrect trials)
% numReps = fr_x(repTrials);                                                    % number of times repeated trials are repeated
% trueNumIncorrs = sum(numReps)-numel(numReps);                                 % check if same as sumReps! Subtract 1 from ea. numReps val for true # incorrects

maxy = max(abs(y));
vel = abs(x);                                                                      % pick up any movement
%    velThold = find(vel > tHold);  NOT IN USE                                     % threshold velocity, lower bound

trialStartInds = find(diff(btTrials) == -1) +1;                                    % transition from 1 to 0; trial start
if excludeTrials == 1
    trialStartInds = trialStartInds(excludeWhichTrials+1:numel(trialStartInds));   % 6.1.16 for initial trial subtractions
    trialStartTimes = time(trialStartInds);
else
    trialStartTimes = time(trialStartInds);
end

% Find proper initial positions (pxy file): will be one interative value before the TS for trialStart
forPosTrialStarts = trialStartInds-1;
posTrialStarts = y(forPosTrialStarts);
%     trueTrials = find(posTrialStarts == posStartsTrueMax |...
%     posStartsTrueMaxNeg | posStartsTrueMin | posStartsTrueMinNeg);          % may need to account for premature movements before last ITI TS

if excludeTrials == 1
    trialEndInds = trialEndInds(excludeTrials+1:numel(trialEndInds));             % 6.1.16 for initial trial subtractions
else
    trialEndInds = find(diff(btTrials) == 1) +1;                                  % these are the last 0s before 1 (before next ITI)
end
trialEndTimes = time(trialEndInds);
PosTrialEnds = y(trialEndInds);
    
%% Assess ALL trial R v L. Note, these are all before long-trial subtractions

trialNum = 1:length(trialStartInds);                                          % absoute number of all trials
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

% sort trials from ITIs - corrected to find the CORRECT bar position at trial start - have to go one i backward for that first value (last ITI)
barpositionCell = {};
m=1; 
% m = excludeTrials+1;
for n=1:2:(numel(trialchange_row)-1);                                           % takes only odd trials (even trials are between trials)
%   barpositioncell{m} = position(trialchange_row(n)+1:trialchange_row(n+1));   % pos data for ea trial put into indvidual cell - Nope, doesn't adjust for mouse moving prematurely
    barpositionCell{m} = position(trialchange_row(n):trialchange_row(n+1));     % takes CORRECTED position data for each trial and puts into indvidual cell
    time_trial{m} = time(trialchange_row(n)+1:trialchange_row(n+1));            % times for each trial - doesn't need to be corrected
    timestart{m} = time(trialchange_row(n)+1);                                  % also doesn't need correcting
    th(m) = position(trialchange_row(n));                                       % added for new "th": bar has to travel abs(th) distance to cross center 
    m=m+1;
end
    
% collect last trial (not detected by diff since no intertrial interval after last trial)
if btTrials(max(trialchange_row)+1) == 0;
    barpositionCell{m} = position(max(trialchange_row):end);                    % takes CORRECTED position data for each trial and puts into indvidual cell
    time_trial{m} = time(max(trialchange_row)+1:end);                           % times for each trial
    timestart{m} = time(max(trialchange_row)+1);
end

%% Subtract initial trials for exclusion:

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

numelLtrials = numel(find(LRtrial<0));
numelRtrials = numel(find(LRtrial>0));

numelRemainingTrials = numelLtrials + numelRtrials;
                     
%% determine which movement is rewarded - updated from Jacki's 2015 version; 2.29.16 - 3.29.16
%   bar starts at neg. position, left WHEEL movement rewarded = -1
%   bar starts at pos. position, right WHEEL movement rewarded = 1

for n = 1:numel(barpositionCell);
    if barpositionCell{n}(1) < 0;                                               % left trials start at (-)bar positions, sorts based on first number of each cell
        reward_side(n)= -1;
    elseif barpositionCell{n}(1) > 0;                                           % right trials start at (+) bar positions, sorts based on first number of each cell
        reward_side(n)= 1;
    end
end

%% set center to define movement as correct "choice" & find incorrects. TS are important here. 
%  correct threshold is distance from starting position to center. BUT if bar hits +/-400 (incorrect "bound"), trial is incorrect
%  BS updated to reflect 2AFC structure. n = 1xnumTrial cell array - corrected for mins (L trial) and maxes (R trial) if animal re-centers bar before solenoid discharges

% Correct L trials:
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) < 0;                                                  % L trial
        if length(find(barpositionCell{n} >= 0) > 0) ;
            thCL(n)= min(find(barpositionCell{n} >= 0));                           % thCL = index of correct threshold crossing for L trials;
            time_thCL(n) = time_trial{n}(min(find(barpositionCell{n} >= 0))-1);    % most importantly, get the TS at thc! (have to go one backward)
        else                                                                       % if bar goes off to left on left trial
            thCL(n) = NaN;  
            time_thCL(n) = NaN;                                                    % NaN is an incorrect L trial; 0 is a R trial or one that was too long
        end                                                                        % if it's finite and > 0, its a correct L trial at that index
    end
end

% what happens to trials that I have to restart for him? - "centered" but no reward
cltNan = find(isnan(time_thCL));
cltZeros = find(time_thCL == 0);
numelCL = numel(time_thCL) - (numel(cltNan) + numel(cltZeros));
percCL = numelCL/numelLtrials * 100;

CLvalid = find(time_thCL > 1);
CLvalidRew = time_thCL(CLvalid);                                                   % centering TS for L correct

% movement extractions for correct L trials
wheelCL = {};
wheelCLtime = {};
wheelLmove = {};
wheelLmoveTime = {};
wheelLfirst = [];
wheelLfirstTime = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) < 0;
        if length(find(barpositionCell{z} >= 0) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) == barpositionCell{z}(i+1) || barpositionCell{z}(i) > barpositionCell{z}(i+1)
                        wheelCL{z}(i) = NaN;
                        wheelCLtime{z}(i) = NaN;
                        wheelLmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelLmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero
                    else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                            wheelCL{z}(i) = barpositionCell{z}(i);              % wheel positions during leftward movements
                            wheelCLtime{z}(i) = time_trial{z}(i);               % all TS of leftward movements
                            wheelLmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST L movement 
                            wheelLmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of L movement
                        end
                    end
                end
            end
        end
    end
end

% initial L bar/wheel movements that lead to centering
for c = 1:length(wheelCL)
    if numel(find(wheelCL{c} >= 0) > 0)
        for i=1:length(wheelCLtime{c})
            if i == length(wheelCLtime{c})
                break
            end
            if wheelLmove{c}(i) > posStartsTrueMaxNeg && wheelLmove{c}(i) < wheelLmove{c}(i+1) && wheelLmove{c}(i) < 0 % still may need to correct for large movements that began < -200
                wheelLfirst(c) = wheelLmove{c}(i);                          % first bar position for L movements 
                wheelLfirstTime(c) = wheelLmoveTime{c}(i);                  % first TS of L movement vectors
                break                                                       % keep i from increasing once it's met this condition
            end
            
        end
    end
end

numelWheelLfirstTime = numel(find(wheelLfirstTime > 1));                    % sanity check: should be same number as numelCL var
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Incorrect L trials. For multiple min values (-400), must find the first one:  updated 10.25.16, 11.12.16
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for n=1:numel(barpositionCell);                                                    % determines position that crosses threshold
    if barpositionCell{n}(1) < 0;                                                  % L trial
        if length(find(barpositionCell{n} <= posStartsTrueMinNegIL) > 0);          % if bar goes off to left on left trial
            for i=1:length(barpositionCell{n})
                find_thIL{n}(i) = barpositionCell{n}(i);
                if length(time_trial{n}) < length(barpositionCell{n})           
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
                if find_thIL{n}(i) <= posStartsTrueMinNegIL
                    thIL(n) = find(find_thIL{n} <= posStartsTrueMinNegIL);
                    time_thIL(n) = find_thILtimes{n}(end);
                    break
                end
                if find_thIL{n}(i) <= posStartsTrueMinNegIL && find_thIL{n}(i+1) <= posStartsTrueMinNegIL
                    find_thIL{n}(i+1) = 0;
                    find_thILtimes{n}(i+1) = 0;
                    thIL(n) = find(find_thIL{n} <= posStartsTrueMinNegIL);
                    time_thIL(n) = find_thILtimes{n}(i(end));
                    break
                end

            end
        else                                                                           % NaN here is a correct L trial; 0 is a R trial
            thIL(n) = NaN;      %this is the first recorded incorrect position (-400)  % if it's finite and > 0, its an incorrect L trial at that index
            time_thIL(n) = NaN; %time of first incorrect position
        end
    end
end

ILvalid = find(time_thIL > 1);     %These are the incorrect L trials         % which trials were incorrect Ls
ILvalidFail = time_thIL(ILvalid);  %And their time stamps!                   % time of bar reaching incorrect threshold

%wheel movement extractions for incorrect L trials
wheelIL = {};
wheelILtime = {};
wheelILmove = {};
wheelILmoveTime = {};
wheelLfirstIL = [];
wheelLfirstTimeIL = [];
clear z;   %for troubleshooting
clear n;

for z=1:numel(barpositionCell)
    if barpositionCell{z}(1) < 0;
        if length(find(barpositionCell{z} <= posStartsTrueMinNegIL)) > 0
            for n = 1:length(find_thIL)                                      % new - correct for only 1 value found for thIL (second condition, below)
                for i=1:length(barpositionCell{z})
                    if barpositionCell{z}(i) > barpositionCell{z}(i+1) && numel(find_thIL{n} > 1)             % moving in incorrect direction
                        wheelIL{z}(i) = barpositionCell{z}(i);               % wheel positions during incorrect L trials
                        wheelILtime{z}(i) = time_trial{z}(i);                % all TS of incorrect L trials
                        wheelILmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST IL movement - Need to change non-position zeroes to NaNs
                        wheelILmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of IL movement
                    end
                    if barpositionCell{z}(i) > barpositionCell{z}(i+1) && numel(find_thIL{n} == 1) % moving in incorrect direction - move this up
                        wheelIL{z}(i) = barpositionCell{z}(i);               % wheel positions during incorrect L trials
                        wheelILtime{z}(i) = time_trial{z}(i);                % all TS of incorrect L trials
                        wheelILmove{z}(i) = barpositionCell{z}(i);           % bit of a hack for the FIRST IL movement - Need to change non-position zeroes to NaNs
                        wheelILmoveTime{z}(i) = time_trial{z}(i);            % & thus, hack for the FIRST TS of IL movement
                        break
                    end
                    if i+1 == length(barpositionCell{z})
                        break
                    else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                            wheelIL{z}(i) = NaN;
                            wheelILtime{z}(i) = NaN;
                            wheelILmove{z}(i) = NaN;                         % prevent repeating positions from being zero
                            wheelILmoveTime{z}(i) = NaN;
                        end
                    end
                end
            end
        end
    end
end

% initial L-trial bar/wheel movements that lead to incorrect threshold
for c = 1:length(wheelILmove)
    for i=1:length(wheelILtime{c})
        if isfinite(wheelILmove{c}(i)) 
            if numel(wheelILmove{c}) == 1
                wheelLfirstIL(c) = wheelILmove{c}(i);                             % first bar position for incorrect L movements
                wheelLfirstTimeIL(c) = wheelILmoveTime{c}(i);                     % first TS of incorrect L movement vectors
             else if numel(wheelILmove{c}) > 1 && wheelILmove{c}(i) < 0
                    wheelLfirstIL(c) = wheelILmove{c}(i);                          % first bar position for incorrect L movements
                    wheelLfirstTimeIL(c) = wheelILmoveTime{c}(i);                  % first TS of incorrect L movement vectors
                    break                                                          % keep i from increasing once it's met this condition
                end
            end
        end
    end
end
numelWheelILfirstTime = numel(find(wheelLfirstTimeIL > 1)); 

%% Correct R trials: 
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) > 0;                                                       % R trial
        if length(find(barpositionCell{n} <= 0) > 0) ;
            thCR(n)= min(find(abs(barpositionCell{n} <= 0)));                           % thCR = index of 1st correct threshold crossing for R trials;
            time_thCR(n) = time_trial{n}(min(find(abs(barpositionCell{n} <= 0)))-1);
        else                                                                            % if bar goes off to left on left trial
            thCR(n) = NaN;              
            time_thCR(n) = NaN;                                                         % NaN is an incorrect R trial; 0 is a L trial
        end                                                                             % so if it's finite and > 0, its a correct R trial at that index
    end
end
crtNan = find(isnan(time_thCR));                                                        % sanity check
crtZeros = find(time_thCR == 0);
numelCR = numel(time_thCR) - (numel(crtNan) + numel(crtZeros));
percCR = numelCR/numelRtrials * 100;

CRvalid = find(time_thCR > 1);
CRvalidRew = time_thCR(CRvalid);                                                        % centering TSs for R correct

% R-trial movement extractions
wheelCR = {};
wheelCRtime = {};
wheelRmove = {};
wheelRmoveTime = {};
wheelRfirst = [];
wheelRfirstTime = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;
        if length(find(barpositionCell{z} <= 0) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) < barpositionCell{z}(i+1) 
                        wheelCR{z}(i) = NaN;
                        wheelCRtime{z}(i) = NaN;
                        wheelRmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelRmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero       
                    else if barpositionCell{z}(i) > barpositionCell{z}(i+1) 
                            wheelCR{z}(i) = barpositionCell{z}(i);              % wheel positions during leftward movements
                            wheelCRtime{z}(i) = time_trial{z}(i);               % all TS of rightward movements
                            wheelRmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST R movement 
                            wheelRmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of R movement
                        end
                    end
                end
            end
        end
    end
end

% initial R-trial movements that lead to centering 
for c = 1:length(wheelCR)
    if numel(find(wheelCR{c} <= 0) > 0)
        for i=1:length(wheelCRtime{c})
            if i == length(wheelCRtime{c})
                break
            end
            if wheelRmove{c}(i) < posStartsTrueMax && wheelRmove{c}(i) > wheelRmove{c}(i+1) && wheelRmove{c}(i) > 0 %still may need to correct for large movements that began > 200
                wheelRfirst(c) = wheelRmove{c}(i);                          % first bar position for R movements
                wheelRfirstTime(c) = wheelRmoveTime{c}(i);                  % first TS of R movement vectors
                break                                                       % keep i from increasing once it's met this condition
            end
            
        end
    end
end
numelWheelRfirstTime = numel(find(wheelRfirstTime > 1));                    % sanity check; should be same as numelCR

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Incorrect R trials:
for n=1:numel(barpositionCell);
    if barpositionCell{n}(1) > 0;                                           % R trial
        if length(find(barpositionCell{n} >= posStartsTrueMaxIR) > 0);      % if bar goes off to right on right trial
            for i=1:length(barpositionCell{n})
                find_thIR{n}(i) = barpositionCell{n}(i);
                if barpositionCell{n}(i) < 0                                % eliminate sneaky correct R trials
                    thIR(n) = NaN;
                    time_thIR(n) = NaN;
                    break
                end
                if length(time_trial{n}) < length(barpositionCell{n})       % adjust for no timestamp even though there is a position at that index
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
                if find_thIR{n}(i) >= posStartsTrueMaxIR
                    thIR(n) = find(find_thIR{n} >= posStartsTrueMaxIR);
                    time_thIR(n) = find_thIRtimes{n}(end);
                    break
                end
                if find_thIR{n}(i) >= posStartsTrueMaxIR && find_thIR{n}(i+1) >= posStartsTrueMaxIR
                    find_thIR{n}(i+1) = 0;
                    find_thIRtimes{n}(i+1) = 0;
                    thIR(n) = find(find_thIR{n} >= posStartsTrueMaxIR);
                    time_thIR(n) = find_thIRtimes{n}(i(end));
                    break
                end
            end
        else                                                                    % NaN here is a correct L trial; 0 is a correct R trial
            thIR(n) = NaN;                                                      % if it's finite and > 0, its an incorrect R trial
            time_thIR(n) = NaN;
        end
    end
end

totCorr = numelCR + numelCL; % num correct trials = 250 (csv); the unaccounted-for Blacktrock trials (missing pulses) = 21; totCorr + 21 = 250

IRvalid = find(time_thIR > 1);
IRvalidFail = time_thIR(IRvalid);

%wheel movement extractions for incorrect R trials
wheelIR = {};
wheelIRtime = {};
wheelIRmove = {};
wheelIRmoveTime = {};
wheelRfirstIR = [];
wheelRfirstTimeIR = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;                                                % So redundant... could totally refine this
        if length(find(barpositionCell{z} >= posStartsTrueMaxIR) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) > barpositionCell{z}(i+1)
                        wheelIR{z}(i) = NaN;
                        wheelIRtime{z}(i) = NaN;
                        wheelIRmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelIRmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero
                  else if barpositionCell{z}(i) < barpositionCell{z}(i+1)        % moving in incorrect direction
                            wheelIR{z}(i) = barpositionCell{z}(i);               % wheel positions during rightward movements
                            wheelIRtime{z}(i) = time_trial{z}(i);                % all TS of incorrect rightward movements
                            wheelIRmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST IR movement 
                            wheelIRmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of IR movement
                        end
                    end
                end
            end
        end
    end
end

% initial R bar/wheel movements that lead to incorrect threshold
clear c;
for c = 1:length(wheelIRmove)
        for i=1:length(wheelIRtime{c})
            if wheelIRmove{c}(i) < posStartsTrueMaxIR && numel(wheelIRmove{c}) < 2 % still may need to correct for large movements that began < -200
                wheelRfirstIR(c) = wheelIRmove{c}(i);                          % first bar position for incorrect R movements 
                wheelRfirstTimeIR(c) = wheelIRmoveTime{c}(i);                  % first TS of incorrect R movement vectors
                break                                                          % keep i from increasing once it's met this condition
            end
            if wheelIRmove{c}(i) < posStartsTrueMaxIR && wheelIRmove{c}(i) < wheelIRmove{c}(i+1) && wheelIRmove{c}(i) > 0  % still may need to correct for large movements that began < -200
                wheelRfirstIR(c) = wheelIRmove{c}(i);                          % first bar position for incorrect R movements 
                wheelRfirstTimeIR(c) = wheelIRmoveTime{c}(i);                  % first TS of incorrect R movement vectors
                break                                                          % keep i from increasing once it's met this condition
            end
        end
end

numelWheelIRfirstTime = numel(find(wheelRfirstTimeIR > 1)); 

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

% centering times:
catCorr = [];
catCorr = horzcat(CLvalidRew, CRvalidRew);
corrects = sort(catCorr');

catRewTrials1 = horzcat(CLvalid, CRvalid);
catRewTrials = sort(catRewTrials1);

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

%Fill a new matrix of rx times:     col 6: CL = -1; CR = 1; IL = 0; IR = 0;
rxMat = NaN(800, 6);
u = 1:totalNumelTrials;
for v = 1:length(u)                 %col1 = trial number
    rxMat(v,1) = u(v);
end
for w = 1:length(rxTimeL)
    if rxTimeL(w) < 0;
        rxMat(w,2) = NaN;
    else
        rxMat(w,2) = rxTimeL(w);    %col2 = rx times for CL trials
        rxMat(w,6) = -1;
    end
end
for x = 1:length(rxTimeR)
    if rxTimeR(x) < 0;
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
        rxMat(y,6) = 0;
    end
end
for z = 1:length(rxTimeIR)
    if rxTimeIR(z) < 0
        rxMat(z,5) = NaN;
    else
        rxMat(z,5) = rxTimeIR(z);    %col2 = rx times for IR trials
        rxMat(z,6) = 0;
    end
end

%bin the rx times:
% hist(x,xbins)
rxCondensed = [];
for m = 1:length(rxMat(:,2))
    if isfinite(rxMat(m,2))
        rxCondensed(m,1) = rxMat(m,1);
        rxCondensed(m,2) = rxMat(m,2);
        rxCondensed(m,3) = rxMat(m,6);
    else if isfinite(rxMat(m,3))
            rxCondensed(m,1) = rxMat(m,1);
            rxCondensed(m,2) = rxMat(m,3);
            rxCondensed(m,3) = rxMat(m,6);
        else if isfinite(rxMat(m,4))
                rxCondensed(m,1) = rxMat(m,1);
                rxCondensed(m,2) = rxMat(m,4);
                rxCondensed(m,3) = rxMat(m,6);
            else if isfinite(rxMat(m,5))
                    rxCondensed(m,1) = rxMat(m,1);
                    rxCondensed(m,2) = rxMat(m,5);
                    rxCondensed(m,3) = rxMat(m,6);
                else rxCondensed(m,1) = rxMat(m,1);
                     rxCondensed(m,2) = NaN;
                     rxCondensed(m,3) = rxMat(m,6);
                end
            end
        end
    end
end

% Sorted Rx times & drag trial info with ea:
rx = rxCondensed(:,2);
trialMat(:,12) = rx(1:length(trialMat)); %to pass Rx info for later trial/move neural alignment; 12th col of trialMat

%% Pull out which wheel L/R firstTimes (correct trials) have fastest/slowest Rx times:

%Fastest Left correct time starts
clear i;
rxCondensed3 = rxCondensed(:,3);
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) < 200
        wheelLfirstTimeFastest(i) = wheelLfirstTime(i);
        end
    end
end
fastestL = wheelLfirstTimeFastest > 1;
numelFastestL = numel(find(fastestL == 1));

%Fastest Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) < 200
        wheelRfirstTimeFastest(j) = wheelRfirstTime(j);
        end
    end
end
fastestR = wheelRfirstTimeFastest > 1;
numelFastestR = numel(find(fastestR == 1));

%mid-rx Left correct time starts
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 200 && rx(i) < 1200;
        wheelLfirstTimeMids(i) = wheelLfirstTime(i);
        end
    end
end

%mid-rx Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 200 && rx(j) < 1200;
        wheelRfirstTimeMids(j) = wheelRfirstTime(j);
        end
    end
end

%slow-rx Left correct time starts
clear i;
for i = 1:length(wheelLfirstTime)
    if i == numel(wheelLfirstTime)
        break
    else if rxCondensed3(i) == -1 && rx(i) >= 1200 && rx(i) < 6000;
        wheelLfirstTimeSlow(i) = wheelLfirstTime(i);
        end
    end
end
if exist('wheelLfirstTimeSlow')
    slowL = wheelLfirstTimeSlow > 1;
    numelSlowL = numel(find(slowL == 1));
else wheelLfirstTimeSlow = [];
end

%mid-rx Right correct time starts
clear j;
for j = 1:length(wheelRfirstTime)
    if j == numel(wheelRfirstTime)
        break
    else if rxCondensed3(j) == 1 && rx(j) >= 1200 && rx(j) < 6000;
        wheelRfirstTimeSlow(j) = wheelRfirstTime(j);
        end
    end
end
if exist('wheelRfirstTimeSlow')
    slowR = wheelRfirstTimeSlow > 1;
    numelSlowR = numel(find(slowR == 1));
else wheelRfirstTimeSlow = [];
end


corrVals = rxCondensed(:,3);
[rxSort1, index] = sort(rx, 1, 'ascend');
rxSort = [];
rxSort(:,1) = index;
rxSort(:,2) = rxSort1;
rxSort(:,3) = corrVals(index);

% edges = [1:5:3000];                  %gives 600 bins for Vgatfive 12.15.16; 23 vals in first bin: numel(find(rxSort(:,2) < 5));
% edges = [5:10:3000];                 %gives 300 bins for Vgatfive 12.15.16; 23 vals
edges = [0:200:5200];                  %gives 27 bins for Vgatfive 12.15.16
% numBins = numel(edges);
mids = conv2(edges, [1 1], 'valid')/2; %gives 26 bins
numBins = numel(edges);

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
mu = rxSort(:,3);
binwidth = edges(2) - edges(1);
binwidthFastest = edges(2) - edges(1);
% corVecFastest = [];
% for i = 1:length(rxSort(:,1))
%     if rxSort1(i) < binwidthFastest 
%         corVecFastest(i) = mu(i);
%     end
% end
% corVecFastest = corVecFastest';
% yesFastest = find(corVecFastest ~= 0); 
% numelCorVecFastest = numel(corVecFastest);
% percFastestCor = numel(yesFastest)/numelCorVecFastest;   %percent corrects for fastest bin
% NumFastest = percFastestCor * 100;
% rxSortFastest = rxSort1(1:binwidthFastest);              %rx times for fastest bin

fastestInds = find(rxSort1 < binwidthFastest);  
rxSortFastest = rxSort1(fastestInds);
% corVecFastest = mu(rxSort1 > edges(1) & rxSort1< binwidthFastest);
corVecFastest = mu(fastestInds);

yesFastest = find(corVecFastest ~= 0); 
numelCorVecFastest = numel(corVecFastest);
percFastestCor = numel(yesFastest)/numelCorVecFastest;   %percent corrects for fastest bin
NumFastest = percFastestCor * 100;

% binwidthFast = mids(4) - mids(1);
binwidthFast = edges(4) - edges(1);
fastInds = find(rxSort1 > edges(2) & rxSort1< binwidthFast);  
rxSortFast = rxSort1(fastInds);                          %rx times for fast bin (200-600 ms)
corVecFast = mu(rxSort1 > edges(1) & rxSort1< binwidthFast);
yesFast = find(corVecFast ~=0);
numelCorVecFast = numel(corVecFast);
percFastCor = numel(yesFast)/numelCorVecFast;
NumFast = percFastCor * 100;                             %percent corrects for fast bin

binwidthMid = edges(6);
midInds = find(rxSort1 > edges(4) & rxSort1 < binwidthMid);  
rxSortMid = rxSort1(midInds);                           %rx times for mid bin (600-1000 ms)
corVecMid = mu(rxSort1 > edges(4) & rxSort1 < binwidthMid);  
yesMid = find(corVecMid ~= 0);
numelCorVecMid = numel(corVecMid);
percMidCor = numel(yesMid)/numelCorVecMid;
NumMid = percMidCor * 100;                              %percent corrects for mid bin

slowInds = find(rxSort1 > binwidthMid);
rxSortSlow = rxSort1(slowInds);                         %rx times for the mid-range bin (600-1000 ms)
corVecSlow = mu(rxSort1 > binwidthMid);
yesSlow = find(corVecSlow ~= 0);
numelCorVecSlow = numel(corVecSlow);
percSlowCor = numel(yesSlow)/numelCorVecSlow;
NumSlow = percSlowCor * 100;                            %percent corrects for slow bin

%Do the plotting
% histRx = histc(rxSort(:,2), edges);                    % all rx times
histRx1 = histc(rxSortFastest, edges);                   % 0-200 ms bin
histRx2 = histc(rxSortFast, edges);                      % 200-600 ms bin
histRx3 = histc(rxSortMid, edges);                       % 600-1000 ms bin
histRx4 = histc(rxSortSlow, edges);                      % >1000 ms bin

figure; hold on;
% h = bar(mids,histRx(1:end-1));                                             
h4 = bar(edges, histRx4, 'k');
h3 = bar(edges, histRx3, 'g'); 
h2 = bar(edges, histRx2, 'b');
h1 = bar(edges, histRx1, 'r'); set(h1, 'BarWidth', 1);  % plot last for visual purposes

    xlabel('Rx times (ms), 200ms bins'); ylabel('# Trials');
%   set(gca, 'Xlim', [min(edges) max(edges)]);
    set(gca, 'Xlim', [0 5200]);
%   set(gca, 'Xtick', [min(edges):600 max(edges)]), %'Xticklabel', [min(edges):max(edges)]);
    set(gca, 'Xtick', [0:500:5200]), 
    num200 = num2str(round(NumFastest)); num200_600 = num2str(round(NumFast)); num600_1000 = num2str(round(NumMid)); num1000 = num2str(round(NumSlow));
    legend(num1000,num600_1000,num200_600,num200,'Location','northeast');
    mouseInfo = strcat('% corrects ', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
    
    
  sumRx1 = vertcat(rxMat(:,2), rxMat(:,3), rxMat(:,4), rxMat(:,5));
 
  meanRx = nanmean(sumRx1);

   
%% Plot the steadiness of correct trials over the session:

% figure(); plot(cumsum(diff(unique([CRvalid CLvalid])))); hold on; plot(diff(cumsum([1:537])));
% percCorrJosh = numel(unique([CRvalid CLvalid]))./537

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
taskbase.wheelLfirstTimeFastest = wheelLfirstTimeFastest;
taskbase.wheelRfirstTimeFastest = wheelRfirstTimeFastest;
taskbase.wheelLfirstTimeMids = wheelLfirstTimeMids;
taskbase.wheelRfirstTimeMids = wheelRfirstTimeMids;
taskbase.wheelLfirstTimeSlow = wheelLfirstTimeSlow;
taskbase.wheelRfirstTimeSlow = wheelRfirstTimeSlow;

taskbase.meanRx = meanRx;


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
    
