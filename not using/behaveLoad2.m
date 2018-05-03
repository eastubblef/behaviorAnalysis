
function [taskbase] = behaveLoad2(fileDirectory)

% 3.1.16 - use as standalone function for structure of task or to be called by recording files for alignment
% USE W/ psychFitTimeTraj.m wrapper (eventually) or TNC_mover.m

% Updated BS 2.28.16 for blinking lum task & 2AFC (modified from bits of psychFitRxnTime3LumsPos3.m and Jacki's older code from July, 2015)
% Output will be only 1 structure with relevant fields "taskbase"
   
%% Initialize

posStartsTrueMax = 200;
posStartsTrueMin = 190;
posStartsTrueMaxNeg = -200; 
posStartsTrueMinNeg = -190;

tHold = 10;     % if velocity is desired - not in use yet
fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks';
filenamestr = 'mVgatthreetag';
    
if exist('filenamestr', 'var');
    [filenamestrE, path] = uigetfile('*.csv*','select the csv file', fileDirectory);                 %BS get the actual # total trials 
end
% 
% if exist('filenamestr', 'var');
%     [filenamestrF, path] = uigetfile('*_p.csv*','select the _p.csv file', fileDirectory);           
% end

if exist('filenamestr', 'var');
    [filenamestrX, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fileDirectory);           %pxy file has solenoid discharge TS with respect to trial starts
end
    
%% Load csv behavior data

clear disp* trial* 
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

%% General information about the file - not quite in use yet 2.29.26; aligns csv and pXY files

%  indArray = 2:size(data.trialTimes,1);                                                  % From the .csv file
%  data.startToEvent  = data.trialTimes(indArray,3) - data.trialTimes(indArray,2);        % use? time from trial start to threshold cross (lever)
%  data.endToStart    = data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4);      % inter-trial interval (should be 3 sec)
%  data.iti = (data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4));              % From the .csv file...I'd rather define it as ITI 
% 
%  data.eventToEnd    = data.trialTimes(indArray,4) - data.trialTimes(indArray,3);        % use? reward delay
%  numTrials   = size(data.trialTimes,1);                                                 % number of total trials
%  lastTrial = data.trialTimes(end,1);                                                    % will be offset by 2 (defined above)

 % BS IMPORTANT: Shift var represents the alignment of continuous data (pXY) with the trial times data (.csv) (eg. state of 0 becomes 2, (last col) & that first col num becomes the trial-start time
 % shift = data.trajectories.contXY(find(data.trajectories.contXY(:,4)>0, 1,'first'),1) - data.trialTimes(1,2); % calculate shift

%      shift = data.trajectories(find(data.trajectories(:,4)>0, 1,'first'),1) - data.trialTimes(1,2); % calculate shift
% 
%      if isempty(shift)
%              shift = 0;
%      end
%          startTimes  = data.trialTimes(:,2) + shift;                        % trial start times
%          eventTimes  = data.trialTimes(:,3) + shift;                        % threshold cross times - relevance for Beth?
%          rewardTimes = data.trialTimes(:,4) + shift;                        % reward times

%% Extract relavent info: use the _pXY file

maxy = max(abs(y));
vel = abs(x);                                                                 % pick up any movement
%    velThold = find(vel > tHold);  NOT IN USE                                % threshold velocity, lower bound

trialStartInds = find(diff(btTrials) == -1) +1;                               % transition from 1 to 0; trial start 
trialStartTimes = time(trialStartInds); 



% Updated 2.28.16 to find proper initial positions (pxy file): will be one interative value before the TS for trialStart
forPosTrialStarts = trialStartInds-1;
posTrialStarts = y(forPosTrialStarts);
%     trueTrials = find(posTrialStarts == posStartsTrueMax |
%     posStartsTrueMaxNeg | posStartsTrueMin | posStartsTrueMinNeg);          % may need to account for premature movements before last ITI TS

trialEndInds = find(diff(btTrials) == 1) +1;                                  % these are the last 0s before 1 (before next ITI)
trialEndTimes = time(trialEndInds);

PosTrialEnds = y(trialEndInds);

%     rewInds = trialEndInds;
%     rewTimes = trialEndTimes;                                               % will be 450-500 ms after centering
    
%% Assess ALL trial R v L 

trialNum = 1:length(trialStartInds);                                          % absoute number of all trials
  lumR = unique(luminance);                                             
  lumL = -unique(luminance);
  lumAll = vertcat(lumR, lumL);
  lumTrials = luminance(trialStartInds);                                      % vector of luminance vals per trial start

%these are ALL of the initial bar/wheel positions
LRtrial = posTrialStarts                                                      % absoute number of total trials
Ltrial = posTrialStarts < 0;
Rtrial = posTrialStarts > 0;

LtrialStarts = trialStartTimes(Ltrial == 1);
LtrialEnds = trialEndTimes(Ltrial == 1);

RtrialStarts = trialStartTimes(Rtrial == 1);
RtrialEnds = trialEndTimes(Rtrial == 1);
    
%% From Jacki's code to get vectors of bar position:     
position = y;                                                                   % bar position (distance for correct centering)
trialChange = diff(btTrials);                                                   % trial changes where zeros change to ones
trialchange_row = find(trialChange);                                            % give row num where trial changes

% sort trials from ITIs - BS corrected to call the CORRECT bar position at trial start - have to go one i backward for that first value (last ITI)
barpositionCell = {};
m=1; 
for n=1:2:(numel(trialchange_row)-1);                                           % takes only odd trials (even trials are between trials)
%   barpositioncell{m} = position(trialchange_row(n)+1:trialchange_row(n+1));   % pos data for ea trial put into indvidual cell - Nope, doesn't adjust for mouse moving prematurely
    barpositionCell{m} = position(trialchange_row(n):trialchange_row(n+1));     % takes CORRECT position data for each trial and puts into indvidual cell
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
    timestart(m) = time(max(trialchange_row)+1);
end

%% determines which movement is rewarded - BS: corrected from Jacki's July 2015 version; updated 2.29.16
%   bar starts at neg. position, left WHEEL movement rewarded = -1
%   bar starts at pos. position, right WHEEL movement rewarded = 1

for n = 1:numel(barpositionCell);
    if barpositionCell{n}(1) < 0;                                               % left trials start at (-)bar positions, sorts based on first number of each cell
        reward_side(n)= -1;
    elseif barpositionCell{n}(1) > 0;                                           % right trials start at (+) bar positions, sorts based on first number of each cell
        reward_side(n)= 1;
    end
end

%% set center to define movement as correct "choice" & find incorrects. TS are importante here. BS updated from Jacki's code 2.29.16 
%  correct threshold is distance from starting position to center. BUT if bar hits +/-400 (incorrect "bound"), trial is incorrect

% BS updated to reflect 2AFC structure. n = 1xnumTrial cell array - corrected for mins (L trial) and maxes (R trial) if animal re-centers bar before solenoid discharges
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

CLvalid = find(time_thCL > 1);
CLvalidRew = time_thCL(CLvalid);                                                   % solenoid TS for L correct

% L wheel movement extractions
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
                else if barpositionCell{z}(i) < barpositionCell{z}(i+1)
                        wheelCL{z}(i) = barpositionCell{z}(i);              % wheel positions during leftward movements
                        wheelCLtime{z}(i) = time_trial{z}(i);               % all TS of leftward movements
                        wheelLmove{z}(i) = barpositionCell{z}(i+1);         % hack for the FIRST L movement
                        wheelLmoveTime{z}(i) = time_trial{z}(i+1);          % & thus, hack for the FIRST TS of L movement
                end       
% %                     if wheelLmove{z}(i) == posStartsTrueMaxNeg || wheelLmove{z}(i) == posStartsTrueMinNeg && wheelLmove{z}(i+1) == 0
% %                             wheelLfirst(z) = find(wheelLmove{z}(i+1));
% %                             wheelLfirstTime(z) = find(time_trial{z}(i+1));
% %                     end    
%                         if wheelLmove{z}(i) > posStartsTrueMaxNeg && wheelLmove{z}(i) < wheelLmove{z}(i+1)
%                            wheelLfirst(z) = wheelLmove{z}(i);        % almost
%                         end 
%                         
%                         if wheelLmoveTime{z}(i) == time_trial{z}(i)
%                             wheelLfirst(z) = wheelLmove{z}(i);
%                         end
% %                     end
                end
            end
        end
    end
end

% wheelCL has positions from posStarts to centering:
% thCL has the indices
thCLnew1 = find(thCL > 0);
thCLnew = thCL(thCLnew1);

% wheelLfirstTime1 = find(wheelLmov

wheelCLtimeValid = [];
wheelCLtimeFind = [];

%% Gotta be a faster/better way:
for c = 1:length(barpositionCell)
    if numel(find(wheelCL{c} >= 0) > 0)
        for i=1:length(wheelCLtime{c})
            for p = 1:length(CLvalidRew)
                if i == length(wheelCLtime{c})
                    break
                end
%                 if i == 1 && wheelCLtime{c}(i) == CLvalidRew(p)             % in case centering occurred on the 1st position, use that time (can't go one back)
%                     wheelCLtimeValid(p) = wheelCLtime{c}(i);
%                     
%                 else if i == 2 && wheelCLtime{c}(i) == CLvalidRew(p)
%                         wheelCLtimeValid(p) = wheelCLtime{c}(i-1);          % in case centering occurred by the 2nd position, take the first
                        
%                     if thCLnew(p) > 2 && wheelCL{c}(i) < -1
% %                              wheelCLtimeFind(p) = min(find(wheelCL{c} > posStartsTrueMaxNeg)); % should find min just > -200
%                             wheelCLtimeFind(p) = min(find(wheelCL{c} < -1));     % should find index of min = -190 or -200
% %                             if min(wheelCL{c}) <= -190 && wheelCL{c}(i) <= -190 && wheelCL{c}(i+1) > -190
%                                if i ==  min(find(wheelCL{c} < -1)) && wheelCLtime{c}(i+1) > wheelCLtime{c}(i)
%                                 wheelCLtimeValid(p) = wheelCLtime{c}(i+1);
%                             end
% 
%                      end
                    
%                     else if wheelCLtime{c}(i) == CLvalidRew(p) && wheelCLtime{c}(i-1) < -1 && wheelCLtime{c}(i-2) == 0
%                             wheelCLtimeValid(p) = wheelCLtime{c}(i-1);
%                         end

                        if wheelLmove{c}(i) > posStartsTrueMaxNeg && wheelLmove{c}(i) < wheelLmove{c}(i+1)
                           wheelLfirst(c) = wheelLmove{c}(i);        % almost
                        end 
                        
                        if wheelLmoveTime{c}(i) == time_trial{c}(i)
                            wheelLfirst(c) = wheelLmove{c}(i);
                        end

                    end
                end
            end
        end
   


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Incorrect L trials. For multiple min values (-400), must correct for finding the first one:
for n=1:numel(barpositionCell);                                                    % determines position that crosses threshold
    if barpositionCell{n}(1) < 0;                                                  % L trial
        if length(find(barpositionCell{n} == -400) > 0);                           % if bar goes off to left on left trial
            for i=1:length(barpositionCell{n})
                find_thIL{n}(i) = barpositionCell{n}(i);
                find_thILtimes{n}(i) = time_trial{n}(i);
                if find_thIL{n}(i) == -400 
                    thIL(n) = find(find_thIL{n} == -400);
                    time_thIL(n) = find_thILtimes{n}(end);
                    break
                end
                
                if find_thIL{n}(i) == -400 && find_thIL{n}(i+1) == -400
                    find_thIL{n}(i+1) = 0;
                    find_thILtimes{n}(i+1) = 0;
                    barpositionCell{n}(i+1) = 0;
                    times_start{n}(i+1) = 0;
                    thIL(n) = find(find_thIL{n} == -400);
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

ILvalid = find(time_thIL > 1);
ILvalidFail = time_thIL(ILvalid);

% Correct R trials: 
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) > 0;                                                       % R trial
        if length(find(barpositionCell{n} <= 0) > 0) ;
            thCR(n)= min(find(abs(barpositionCell{n} <= 0)));                           % thCR = index of 1st correct threshold crossing for R trials;
            time_thCR(n) = time_trial{n}(min(find(abs(barpositionCell{n} <= 0)))-1);    % most importantly, get the TS at thc! (have to go one backward)
        else                                                                            % if bar goes off to left on left trial
            thCR(n) = NaN;              
            time_thCR(n) = NaN;                                                         % NaN is an incorrect R trial; 0 is a L trial
        end                                                                             % so if it's finite and > 0, its a correct R trial at that index
    end
end
crtNan = find(isnan(time_thCR));                                                        % sanity check
crtZeros = find(time_thCR == 0);
numelCR = numel(time_thCR) - (numel(crtNan) + numel(crtZeros));

CRvalid = find(time_thCR > 1);
CRvalidRew = time_thCR(CRvalid);                                                        % solenoid TS for R correct

% wheel movement extractions to the right
for z=1:numel(barpositionCell);
    if barpositionCell{z}(1) > 0;                                                  
        if length(find(barpositionCell{z} >= 0) > 0)
            for i=1:length(barpositionCell{z})
                if i+1 == length(barpositionCell{z})
                    break
                else if barpositionCell{z}(i) > barpositionCell{z}(i+1)
                        wheelCR{z}(i) = barpositionCell{z}(i);
                        wheelCRtime{z}(i) = time_trial{z}(i);
                        
                    end
                end
            end
        end
    end
end

% Incorrect R trials:
for n=1:numel(barpositionCell);                                                 
    if barpositionCell{n}(1) > 0;                                                       % R trial
        if length(find(barpositionCell{n} == 400) > 0);                                 % if bar goes off to left on left trial
            for i=1:length(barpositionCell{n})
                find_thIR{n}(i) = barpositionCell{n}(i);                                % correct for multiple maxes == 400
                find_thIRtimes{n}(i) = time_trial{n}(i);
                if find_thIR{n}(i) == 400 
                    thIR(n) = find(find_thIR{n} == 400);
                    time_thIR(n) = find_thIRtimes{n}(end);
                    break
                end
                if find_thIR{n}(i) == 400 && find_thIR{n}(i+1) == 400
                   find_thIR{n}(i+1) = 0;
                    find_thIRtimes{n}(i+1) = 0;
                    barpositionCell{n}(i+1) = 0;
                    times_start{n}(i+1) = 0;
                    thIR(n) = find(find_thIR{n} == -400);
                    time_thIR(n) = find_thIRtimes{n}(i(end));
                    break
                end
            end
        else                                                                        % NaN here is a correct L trial; 0 is a correct R trial
            thIR(n) = NaN;                                                          % if it's finite and > 0, its an incorrect R trial
            time_thIR(n) = NaN;
        end
    end
end

totCorr = numelCR + numelCL; %AHA! total num correct trials = 250 (csv file); the unaccounted-for Blacktrock trials (missing pulses) = 21; totCorr + 21 = 250! YES!!

IRvalid = find(time_thIR > 1);
IRvalidFail = time_thIR(IRvalid);

%% Tally up any trials that started & never were rewarded or incorrect (due to mouse movement and wheel centering b4 trial start)
csvTrials = data.attemptedTrials;
pxyTrials = trialStartTimes;
diffTrials = numel(pxyTrials) - numel(csvTrials);

%Would be nice to get a matrix of trialStartTimes (Lcol) and trialEvents(Rcol) NaN = not completed trial = number of diffTrials
catCorr = [];
catCorr = horzcat(CLvalidRew, CRvalidRew);
corrects = sort(catCorr');

% trialStartTimes = trialStartTimes';
TSstartsNcorrects = vertcat(trialStartTimes, corrects);
sortTSstartsNcorrects = sort(TSstartsNcorrects);

% need solenoid clicks:
rewTimes = [];
for t = 1:length(trialEndTimes);
    for c = 1:length(corrects);                                                % solenoid click will happen bt t and t+1 for trial t & first trial is always rewarded
        if corrects(c) > trialEndTimes(t) && corrects(c) < trialEndTimes(t+1) 
            rewTimes(c-1) = trialEndTimes(t);
        end
    end
end
  
startsRewCat = horzcat(trialStartTimes', rewTimes);
startsRews = sort(startsRewCat);
% validStrtsRews = find(diff(startsRews)>50);  % account for very low vals: rx time that is virtually impossible due to mouse likely already moving
% validStrtsRews = startsRews(validStrtsRews);

%% Create taskbase structure (rows)

taskbase.rewTimes = rewTimes;
taskbase.csvTrialInds = csvTrials';
taskbase.LtrialStarts = LtrialStarts;
taskbase.RtrialStarts = RtrialStarts;
taskbase.leftCenteredTimes = CLvalidRew;
taskbase.leftIndsCenteredTimes = time_thCL;
taskbase.leftIncorrectTimes = ILvalidFail;
taskbase.rightCenteredTimes = CRvalidRew;
taskbase.rightIndsCenteredTimes = time_thCR;
taskbase.rightIncorrectTimes = IRvalidFail;
taskbase.trialStartInds = trialStartInds';
taskbase.trialStartTimes = trialStartTimes';
taskbase.diffTrials = diffTrials;                                               % if below number doesn't = num Blackrock pulses, add this num to it
taskbase.events = sortTSstartsNcorrects';                                       % numel here should = num pulses into blackrock for .nev file
taskbase.startsRews = startsRews;
taskbase.LRtrial = LRtrial;

%% COMPLETED ALL ANALYSIS. SAVE taskase STRUCTURE

disp(' ');    
disp('%-------------------------------------------------------------------');
disp(['Completed file: ' filenamestrX]);
  
filenamestrNew = fileDirectory;
targetName = filenamestrX(1:end-11);     
cd(filenamestrNew);
     
 save([targetName '_tb.mat'],'taskbase');                              
 disp(['saved as ' targetName '_tb.mat']);
                                                                
%     save(FR_file{n_i}, 'FR', 'FR_bine');
%     save([targetName '_bh.mat'],'ContData');      
     
disp('%-------------------------------------------------------------------');
disp(' ');
    
