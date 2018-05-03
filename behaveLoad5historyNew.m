
function [taskbase] = behaveLoad5historyNew(fileDirectory)
% 4.27.16 - updated for more reliable incorrect L/R trial extraction - THIS IS THE LATEST VERSION
% 4.16.16 - updated lines 397 & 398 to accomodate behavioral session 151106
% 3.28.16 - updated for initial incorrect L and R movements 
% 3.24.16 - updated for initial L and R movements that lead to reward
% 3.1.16 - use as standalone function for structure of task or to be called by recording files for alignment
% USE W/ psychFitTimeTraj.m wrapper (eventually) or TNC_mover.m

% Updated BS 2.28.16 for blinking lum task & 2AFC (modified from bits of psychFitRxnTime3LumsPos3.m and Jacki's older code from July, 2015)
% Output will be only 1 structure with relevant fields "taskbase"

%8.15.16:

   
%% Initialize
posStartsTrueMaxIR = 400;
% posStartsTrueMaxIR = 510;
posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMinNeg = -190;
posStartsTrueMaxNeg = -200; 
posStartsTrueMinNegIL = -400;
% posStartsTrueMinNegIL = -510;

% %blocked trial constraints:
% trialConstrainMin = 60;
% trialConstrainMax = 170; %best performance within 100ish trials


tHold = 10;     % if velocity is desired - not in use yet
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151118/behaveChunks';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgattwo/151119/151119good/behaveChunks/untagged';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151106/behaveChunks';
% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151104/1001/behaveChunks';
% fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/Vgatfive/forAnalysis';
 fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/Vgatfive/forAnalysis';

% filenamestr = 'mVgatthree';
% filenamestr = 'mVgattwo'; 
filenamestr = 'mVgatfive';
% filenamestr = 'mVgatsix';

if exist('filenamestr', 'var');
    [filenamestrE, path] = uigetfile('*.csv*','select the csv file', fileDirectory);                 % get the actual # total trials 
end

% if exist('filenamestr', 'var');
%     [filenamestrF, path] = uigetfile('*_p.csv*','select the _p.csv file', fileDirectory);           
% end

if exist('filenamestr', 'var');
    [filenamestrX, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fileDirectory);          % pxy file has solenoid discharge TS with respect to trial starts
end
    
%% Load csv behavior data

% clear disp* trial* 
 skipInitTrials = 0;

 csvFile = filenamestrE;
 data.trialTimes = dlmread(csvFile,',',2,0);

 data.attemptedTrials = data.trialTimes(:,1);
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

maxy = max(abs(y));
vel = abs(x);                                                                 % pick up any movement
%    velThold = find(vel > tHold);  NOT IN USE                                % threshold velocity, lower bound

trialStartInds = find(diff(btTrials) == -1) +1;                               % transition from 1 to 0; trial start 
trialStartTimes = time(trialStartInds); 

% Updated 2.28.16 to find proper initial positions (pxy file): will be one interative value before the TS for trialStart
forPosTrialStarts = trialStartInds-1;
posTrialStarts = y(forPosTrialStarts);
%     trueTrials = find(posTrialStarts == posStartsTrueMax |...
%     posStartsTrueMaxNeg | posStartsTrueMin | posStartsTrueMinNeg);          % may need to account for premature movements before last ITI TS

trialEndInds = find(diff(btTrials) == 1) +1;                                  % these are the last 0s before 1 (before next ITI)
trialEndTimes = time(trialEndInds);
PosTrialEnds = y(trialEndInds);
    
%% Assess ALL trial R v L 

trialNum = 1:length(trialStartInds);                                          % absoute number of all trials
lumR = unique(luminance);                                             
lumL = -unique(luminance);
lumAll = vertcat(lumR, lumL);
lumTrials = luminance(trialStartInds);                                        % vector of luminance vals per trial start

%these are ALL of the initial bar/wheel starting positions
LRtrial = posTrialStarts;                                                     % absoute number of total trials
Ltrial = posTrialStarts < 0;
Rtrial = posTrialStarts > 0;

numelLtrials = numel(find(LRtrial<0));
numelRtrials = numel(find(LRtrial>0));

LtrialStarts = trialStartTimes(Ltrial == 1);
if numel(Ltrial) > numel(trialEndTimes)
    Ltrial = Ltrial(end-1);
end
LtrialEnds = trialEndTimes(Ltrial == 1);

RtrialStarts = trialStartTimes(Rtrial == 1);
if numel(Rtrial) > numel(trialEndTimes)
    Rtrial = Rtrial(end-1);
end
RtrialEnds = trialEndTimes(Rtrial == 1);
    
%% From Jacki's code to get vectors of bar position:     
position = y;                                                                   % bar position (distance for correct centering)
trialChange = diff(btTrials);                                                   % trial changes where zeros change to ones
trialchange_row = find(trialChange);                                            % give row num where trial changes

% sort trials from ITIs - corrected to find the CORRECT bar position at trial start - have to go one i backward for that first value (last ITI)
barpositionCell = {};
m=1; 
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
%     barpositioncell{m} = position(max(trialchange_row)+1:end);                % pos data for ea trial put into cell - Nope, doesn't adjust for mouse moving prematurely
    barpositionCell{m} = position(max(trialchange_row):end);                    % takes CORRECTED position data for each trial and puts into indvidual cell
    time_trial{m} = time(max(trialchange_row)+1:end);                           % times for each trial    
    timestart{m} = time(max(trialchange_row)+1);
end

%% determines which movement is rewarded - updated from Jacki's 2015 version; 2.29.16 - 3.29.16
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
            time_thCL(n) = NaN;                                                    % NaN is an incorrect L trial; 0 is a R trial
        end                                                                        % if it's finite and > 0, its a correct L trial at that index
    end
end
% numelCL = numel(find(isfinite(time_thCL > 0)));                                  % sanity check
cltNan = find(isnan(time_thCL));
cltZeros = find(time_thCL == 0);
numelCL = numel(time_thCL) - (numel(cltNan) + numel(cltZeros));
percCL = numelCL/numelLtrials * 100;

CLvalid = find(time_thCL > 1);
CLvalidRew = time_thCL(CLvalid);                                                   % centering TS for L correct

% L wheel movement extractions for correct L trials
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
% Incorrect L trials. For multiple min values (-400), must find the first one:  updated 4.27.16
for n=1:numel(barpositionCell);                                                    % determines position that crosses threshold
    if barpositionCell{n}(1) < 0;                                                  % L trial
        if length(find(barpositionCell{n} <= posStartsTrueMinNegIL) > 0);          % if bar goes off to left on left trial
            for i=1:length(barpositionCell{n})
                find_thIL{n}(i) = barpositionCell{n}(i);
                
                 if length(time_trial{n}) < length(barpositionCell{n})
                    time_trial{n}(i) = NaN;
                 end

                
%                 find_thILtimes{n}(i) = time_trial{n}(i);
                
                if isempty(time_trial{n}(i))                                       % in case timestamp wasn't collected (160516)
                    find_thILtimes{n}(i) = time_trial{n}(i-1);
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
        else                                                                        % NaN here is a correct L trial; 0 is a R trial
            thIL(n) = NaN;                                                          % if it's finite and > 0, its an incorrect L trial at that index
            time_thIL(n) = NaN;
        end
    end
end

ILvalid = find(time_thIL > 1);                                                      % which trials were incorrect Ls
ILvalidFail = time_thIL(ILvalid);                                                   % time of bar reaching incorrect threshold

% L wheel movement extractions for incorrect L trials
wheelIL = {};
wheelILtime = {};
wheelILmove = {};
wheelILmoveTime = {};
wheelLfirstIL = [];
wheelLfirstTimeIL = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) < 0;
        if length(find(barpositionCell{z} <= posStartsTrueMinNegIL) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                        wheelIL{z}(i) = NaN;
                        wheelILtime{z}(i) = NaN;
                        wheelILmove{z}(i) = NaN;                                 % prevent repeating positions from being zero
                        wheelILmoveTime{z}(i) = NaN;                             % prevent repeating time(positions) from being zero
                  else if barpositionCell{z}(i) > barpositionCell{z}(i+1)        % moving in incorrect direction
                            wheelIL{z}(i) = barpositionCell{z}(i);               % wheel positions during leftward movements
                            wheelILtime{z}(i) = time_trial{z}(i);                % all TS of incorrect leftward movements
                            wheelILmove{z}(i) = barpositionCell{z}(i+1);         % bit of a hack for the FIRST IL movement - Need to change non-position zeroes to NaNs
                            wheelILmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of IL movement
                        end
                    end
                end
            end
        end
    end
end

% initial L bar/wheel movements that lead to incorrect threshold
for c = 1:length(wheelILmove)
    if numel(find(wheelILmove{c} <= posStartsTrueMinNegIL) > 0)
        for i=1:length(wheelILtime{c})
            if i == length(wheelILtime{c})
                break
            else if wheelILmove{c}(i) < posStartsTrueMaxNeg && wheelILmove{c}(i) > wheelILmove{c}(i+1) && wheelILmove{c}(i) < 0 % still may need to correct for large movements that began < -200
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
%             time_thCR(n) = time_trial{n}(min(find(abs(barpositionCell{n} <= 0)))-1);    % most importantly, get the TS at thc! (have to go one backward)
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

% R wheel movement extractions
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

% initial R bar/wheel movements that lead to centering 
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
                % if i > length(barpositionCell{n})
                if i > length(time_trial{n})  %added for 151106
                    break
                else if length(time_trial{n} > 1)
                        find_thIR{n}(i) = barpositionCell{n}(i);            % correct for multiple maxes == 400
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
            end
        else                                                                % NaN here is a correct L trial; 0 is a correct R trial
            thIR(n) = NaN;                                                  % if it's finite and > 0, its an incorrect R trial
            time_thIR(n) = NaN;
        end
    end
end

totCorr = numelCR + numelCL; % num correct trials = 250 (csv); the unaccounted-for Blacktrock trials (missing pulses) = 21; totCorr + 21 = 250

IRvalid = find(time_thIR > 1);
IRvalidFail = time_thIR(IRvalid);

% R wheel movement extractions for incorrect R trials
wheelIR = {};
wheelIRtime = {};
wheelIRmove = {};
wheelIRmoveTime = {};
wheelRfirstIR = [];
wheelRfirstTimeIR = [];
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;
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
for c = 1:length(wheelIRmove)
    if numel(find(wheelIRmove{c} >= posStartsTrueMaxIR) > 0)
        for i=1:length(wheelIRtime{c})
            if i == length(wheelIRtime{c})
                break
            end
            if wheelIRmove{c}(i) < posStartsTrueMaxIR && wheelIRmove{c}(i) < wheelIRmove{c}(i+1) && wheelIRmove{c}(i) > 0 % still may need to correct for large movements that began < -200
                wheelRfirstIR(c) = wheelIRmove{c}(i);                          % first bar position for incorrect R movements 
                wheelRfirstTimeIR(c) = wheelIRmoveTime{c}(i);                  % first TS of incorrect R movement vectors
                break                                                          % keep i from increasing once it's met this condition
            end
        end
    end
end

numelWheelIRfirstTime = numel(find(wheelRfirstTimeIR > 1));                    

%% Tally up any trials that started & never were rewarded or incorrect (due to mouse movement and wheel centering b4 trial start)
csvTrials = data.attemptedTrials;
pxyTrials = trialStartTimes;
diffTrials = numel(pxyTrials) - numel(csvTrials);

% centering times
catCorr = [];
catCorr = horzcat(CLvalidRew, CRvalidRew);
corrects = sort(catCorr');

TSstartsNcorrects = vertcat(trialStartTimes, corrects);
sortTSstartsNcorrects = sort(TSstartsNcorrects);

% solenoid clicks
rewTimes = [];
for t = 1:length(trialEndTimes);
    for c = 1:length(corrects);                                             % solenoid click will happen bt t and t+1 for trial t; first trial is always rewarded
        if corrects(c) > trialEndTimes(t) && corrects(c) < trialEndTimes(t+1) 
%             rewTimes(c-1) = trialEndTimes(t);
            rewTimes(c) = trialEndTimes(t+1);    % for 160505 Vgatfour dataset

        end
    end
end
  
percCorrectsTot = totCorr/numel(lumTrials) * 100;


startsRewCat = horzcat(trialStartTimes', rewTimes);
startsRews = sort(startsRewCat);
% validStrtsRews = find(diff(startsRews)>50);                               % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
% validStrtsRews = startsRews(validStrtsRews);

% sanity matrix (correct trials), 6 cols: trialStarts, Lfirst rewarded movement, Rfirst rewarded movement, Lcentered, Rcentered, trial end times
% rxTimeL = [];
% rxTimeR = [];
% wheelLfirstTime = wheelLfirstTime';
% wheelRfirstTime = wheelRfirstTime';
% trialMat = NaN(length(trialEndTimes), 6);
% if numel(trialEndTimes) < numel(trialStartTimes)
%     trialStartTimesNew = trialStartTimes(end-1);
% end
% 
% trialMat(:,1) = trialStartTimes;
% trialMat(:,6) = trialEndTimes;
% for i = 1:length(wheelLfirstTime)
%         trialMat(i,2) = wheelLfirstTime(i);
%         rxTimeL(i) = wheelLfirstTime(i) - trialStartTimes(i);
%         trialMat(i,4) = time_thCL(i);
% end
% for i = 1:length(wheelRfirstTime)
%         trialMat(i,3) = wheelRfirstTime(i);
%         rxTimeR(i) = wheelRfirstTime(i) - trialStartTimes(i);
%         trialMat(i,5) = time_thCR(i);
% end
%  
% rxL = rxTimeL > 1;
% rxLvalid = rxTimeL(rxL);
% meanRxLvalid = mean(rxLvalid);
% medianRxLvalid = median(rxLvalid);
% 
% rxR = rxTimeR > 1;
% rxRvalid = rxTimeR(rxR);
% meanRxRvalid = mean(rxRvalid);
% medianRxRvalid = median(rxRvalid);

%% Create ability to look at individual trial history 8.10.16
% for all L choices, concatenate CL and IR trials:

totLchoiceTrials = sort(horzcat(CLvalid, IRvalid));
percLchoiceTot = numel(totLchoiceTrials)/n * 100;

%plot actual trial structure:
LRtrialStruct = [];
for i = 1:length(LRtrial)
    if LRtrial(i) < 1
        LRtrialStruct(i) = -1;
    else if LRtrial(i) > 1
            LRtrialStruct(i) = 1;
        end 
    end
end

%plot actual choice:
choiceVec = [];
for i = 1:length(trialNum)
    if any(find(trialNum(i) == totLchoiceTrials))
        choiceVec(i) = -1;
    else
        choiceVec(i) = 1;
    end
end


% %% For L choice, 1 trial back, 2 trials back... up to i-7 (trials 60-170 are the constraints set):
% 
% choiceVecConstrained = choiceVec(trialConstrainMin:trialConstrainMax);          %meaningful vector of choices from constrained start to finish
% 
% %Define trial i:
% Li = numel(choiceVecConstrained);
% prob_L = find(choiceVecConstrained < 0);                %these are all the L choices within the constrained trials
% numProbL = numel(prob_L);
% Prob_Li = numProbL/Li * 100;                            %This will be the first x value on P(Li) at i
% 
% %Define trial i-1:
% Lis = [];
% for i = 2:length(choiceVecConstrained)
%     Lis(i) = choiceVecConstrained(i-1) < 0;              %find the -1s occuring at i-1
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% % prob_LminusOne = choiceVecConstrained > 0;
% prob_LminusOne = Lis > 0;
% numProbLminusOne = find(prob_LminusOne == 1);
% Prob_LiMinusOne = numel(numProbLminusOne)/(Li-1) * 100;  %Must divide by true numel (i=2:end); %This will be the second x value on P(Li) at i-1
% 
% %Define trial i-2:
% Lis2 = [];
% for i = 3:length(choiceVecConstrained)
%     Lis2(i) = choiceVecConstrained(i-2) < 0;             %find the -1s occuring at i-2
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusTwo = Lis2 > 0;
% numProbLminusTwo = find(prob_LminusTwo == 1);
% Prob_LiMinusTwo = numel(numProbLminusTwo)/(Li-2) * 100;
% 
% %Define trial i-3:
% Lis3 = [];
% for i = 4:length(choiceVecConstrained)
%     Lis3(i) = choiceVecConstrained(i-3) < 0;             %find the -1s occuring at i-3
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusThree = Lis3 > 0;
% numProbLminusThree = find(prob_LminusThree == 1);
% Prob_LiMinusThree = numel(numProbLminusThree)/(Li-3) * 100;
% 
% %Define trial i-4:
% Lis4 = [];
% for i = 5:length(choiceVecConstrained)
%     Lis4(i) = choiceVecConstrained(i-4) < 0;             %find the -1s occuring at i-3
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusFour = Lis4 > 0;
% numProbLminusFour = find(prob_LminusFour == 1);
% Prob_LiMinusFour = numel(numProbLminusFour)/(Li-4) * 100;
% 
% %Define trial i-5:
% Lis5 = [];
% for i = 6:length(choiceVecConstrained)
%     Lis5(i) = choiceVecConstrained(i-5) < 0;             %find the -1s occuring at i-3
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusFive = Lis5 > 0;
% numProbLminusFive = find(prob_LminusFive == 1);
% Prob_LiMinusFive = numel(numProbLminusFive)/(Li-5) * 100;
% 
% %Define trial i-6:
% Lis6 = [];
% for i = 7:length(choiceVecConstrained)
%     Lis6(i) = choiceVecConstrained(i-6) < 0;             %find the -1s occuring at i-3
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusSix = Lis6 > 0;
% numProbLminusSix = find(prob_LminusSix == 1);
% Prob_LiMinusSix = numel(numProbLminusSix)/(Li-6) * 100;
% 
% %Define trial i-7:
% Lis7 = [];
% for i = 8:length(choiceVecConstrained)
%     Lis7(i) = choiceVecConstrained(i-7) < 0;             %find the -1s occuring at i-3
%     if i == choiceVecConstrained(end)
%         break
%     end
% end
% prob_LminusSeven = Lis7 > 0;
% numProbLminusSeven = find(prob_LminusSeven == 1);
% Prob_LiMinusSeven = numel(numProbLminusSeven)/(Li-7) * 100;

%% Plot Prob(left trial i) as a function of trial history:
% figure; 
%     hold on;     
% %     gcf; 
% % hold on; 
% y = [Prob_Li_jitter, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1];        
% h = plot(y, 'r.', 'MarkerSize', 25, 'Color', [0.5 0.5 0.5]);
%     ylabel('P choice (Left_i)', 'FontSize', 22);  
%     xlabel('Trial history (i-num trials back)', 'FontSize', 22); 
%     ylim([0 1.05]);
%     set(gca, 'YTick', 0:.20:1.2, 'XTick', 0:1:8);
%     ax = gca; 
%     set(ax);
%     ax.XTickLabelMode = 'manual'
%     ax.XTickLabel = {' ', 'i','i-1','i-2','i-3', 'i-4','i-5','i-6','i-7'}
%     ax.FontSize = 22;


%% Create taskbase structure (row vectors)

taskbase.rewTimes = rewTimes;
taskbase.csvTrialInds = csvTrials';
taskbase.LtrialStarts = LtrialStarts';
taskbase.RtrialStarts = RtrialStarts';
taskbase.leftCenteredTimes = CLvalidRew;
taskbase.leftIndsCenteredTimes = time_thCL;
taskbase.leftIncorrectTimes = ILvalidFail;
taskbase.rightCenteredTimes = CRvalidRew;
taskbase.rightIndsCenteredTimes = time_thCR;
taskbase.rightIncorrectTimes = IRvalidFail;
taskbase.trialStartInds = trialStartInds';
taskbase.trialStartTimes = trialStartTimes';
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

taskbase.percCorrectsTot = percCorrectsTot;
taskbase.choiceVec = choiceVec;
taskbase.LRtrialStruct = LRtrialStruct;

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
    
