function [taskbase_bar] = excel2taskbase_bar (rawdatafile)
 %% Updated 2.28.16 BS for 2AFC task, bound of 200, and trial-aligment with recording clock

%     if nargin < 3, bound = 200; % when bound unspecified, set to 200      % no longer necessary since this is in the pxy file
%     end                                                                   % see position var below:

%   excel2taskbase_bar_v1 - transforms data format from excel (.csv) to
%   taskbase_bar structure
 
%   Reads data in files saved from excel(.csv files) and transforms it 
%   to a user-friendly 'taskbase_bar' data format.
%   The 'taskbase_bar' format is meant to mirror that created by disp2taskbase.
% 
%   Currently works for bar task data saved in excel format
%
%   This code extracts from pXY.csv files
%   
%   INPUT
%   rawdatafile = .csv file 
%   th = threshold
%   bound = maximum eccentricity of the bar 
%   
%   OUTPUT
%   taskbase_bar:
%
%   taskbase_bar.trajectory =             barpositioncell;  all trial trajectories (cell array)
%   taskbase_bar.trial_start =            timestart; start times for all trials(Nx1 vector)
%   taskbase_bar.luminance =              brightnesscell; lum of bar for all trials (Nx1 vector)
%   taskbase_bar.reward_side =            reward_side; left=-1, right=1 (Nx1 vector)
%   taskbase_bar.choice =                 choice; actual side chosen; left=-1, right=1, no choice = 0 (Nx1 vector)
%                                           update for 2AFC task
%   taskbase_bar.initiation =             initiation; time of 'choice' movement start(Nx1 vector)
%   taskbase_bar.expid =                  rawdatafile; 
%   expid =                                experiment id; gives file name
%   taskbase_bar.displacement_threshold = th; absolute movement from starting position needed to be considered a choice (1-element vector)
%   taskbase_bar.displacement_bound =     bound; absolute limit position that the bar can move
%                                           away from target (1-element vector)

%% define matrix with bar position, time, bar lum and bar direction as "bartaskdata"
%   reads pXY.csv file starting at the second third row and first column
%   files: Vgatone_2015_03_17_1627pXY.csv;
%          Gad2thr_2014_12_05_1442_pXY.csv;
%          Vgatone_2015_03_18_178pXY.csv;

%   BS modified 2.28.16
%   BS modified 7/3/15

%% 
bartaskdata = csvread(rawdatafile,2,0); %starts reading at row 2, column zero (this is where the values start)
%   index each variable in matrix

position = bartaskdata(:,3);    % initial bar position and how far bar must go for centering
% th = position;                  % BS updated for 2AFC task - th is center; thus th is position as a distance
timeMS = bartaskdata(:,1);      % preserves ms timesscale
time =(bartaskdata(:,1)/1000);  % changes time from msec to sec
trial = bartaskdata(:,5);       % each trial is either all 0's; intertrial is all 1's
brightness = bartaskdata(:,7);  % BS updated for blinking bar 
trialchange = diff(trial);      % determine where trial changes based on when 0's change to 1's and vice versa
trialchange_row = find(trialchange); % gives row number where trial changes

%% sort trials from intertrial intervals - BS: Be aware, this doesn't call the CORRECT bar position at trial start
barpositioncell = {}; % empty cell array for bar position data

m=1; 
for n=1:2:(numel(trialchange_row)-1);                                         % takes only odd trials (even trials are between trials)
    barpositioncell{m} = position(trialchange_row(n)+1:trialchange_row(n+1)); % pos data for ea trial put into indvidual cell - Nope, doesn't adjust for mouse moving prematurely
%     th{m} = barpositioncell{1,n}(1);
    time_trial{m} = time(trialchange_row(n)+1:trialchange_row(n+1));          % times for each trial
    timestart(m) = time(trialchange_row(n)+1);
    th(m) = position(trialchange_row(n)+1);     %must update th to reflect true position of bar BEFORE trial start (can only be 200 or 190)
    m=m+1;
end

% th = barpositioncell{1,1}(1); %use this to index into the first value of each barposition cell

%   collects last trial (not detected by diff since no intertrial interval after last trial)
if trial(max(trialchange_row)+1) == 0;
    barpositioncell{m} = position(max(trialchange_row)+1:end); % takes position data for each trial and puts into indvidual cell
    time_trial{m} = time(max(trialchange_row)+1:end);          % times for each trial
    timestart(m) = time(max(trialchange_row)+1);
end



%% lum for each trial
%   gives luminance of bar for each trial & th
m=1;
for n=1:2:numel(trialchange_row);
    brightness_trial(m) = brightness(trialchange_row(n));                   % gives corresponding bar brightness for each trial
    th(m) = barposition
    m=m+1;
end

%% determines which movement is rewarded - BS: incorrect. 2.29.16
%   bar starts at pos. position, left WHEEL movement rewarded= -1
%   bar starts at neg. position, right WHEEL movement rewarded= 1

for n = 1:numel(barpositioncell);
    if barpositioncell{n}(1) < 0;     % left trials start at (-)bar positions, sorts based on first number of each cell
        reward_side(n)= 1;
    elseif barpositioncell{n}(1) > 0; % right trials start at (+) bar positions, sorts based on first number of each cell
        reward_side(n)= -1;
    end
end

%% set threshold to define movement as a "choice"
%   threshold is absolute distance from starting position
%   HOWEVER, if the bound is hit (|bound|), that is considered a choice
%
%      left choice = -1
%      right choice = 1
%      no choice = 0

%   finds if threshold is crossed
%   th = threshold; bar needs to be moved at least threshold distance from start to be a choice; part of function input
%   bound = absolute limit that the bar can move away from target; part of function input 

% BS : work on this to reflect 2AFC structure. n = 1xnumTrial cell array
for n=1:numel(barpositioncell);                                                % determines position that crosses threshold
    if length(find(abs(barpositioncell{n}-barpositioncell{n}(1))> th)) > 0;
        thc(n)= min(find(abs(barpositioncell{n}-barpositioncell{n}(1)) > th)); % thc=threshold crossing;
    else
        thc(n) = 0;                                                            % if threshold isn't crossed, "position" is zero (could also be NaN)
    end
end

%   determine choice (left, right, or no choice) based on where threshold is crossed
m=1;
for n=1:numel(barpositioncell);
    if thc(n)== 0;
        choice(n) = 0;% no choice; no movements met threshold
    elseif any(barpositioncell{n}(1): barpositioncell{n}(thc(n)) == -1*position & barpositioncell{n}(1)~= -1*bound); % neg. bound hit; WHEEL moved left
        choice(n) = -1;
    elseif barpositioncell{n}(thc(n)) < barpositioncell{n}(1); % end position less than starting position; WHEEL moved left
        choice (n) = -1;
    elseif any(barpositioncell{n}(1): barpositioncell{n}(thc(n)) == bound & barpositioncell{n}(1) ~= bound);% pos. bound hit; WHEEL moved right
        choice(n) = 1;
    elseif barpositioncell{n}(thc(n)) > barpositioncell{n}(1);% end position greater than starting position; WHEEL moved right
        choice(n) = 1;
    end
end



%% find threshold movement initition time
%   finds time of first movement that is in the same direction of the threshold crossing movement
%   times are running time of the session

for n = 1:numel(choice);
    if choice(n) == 0;
        initiation(n) = NaN;                   % considers time of first movement as zero
    elseif choice(n) == -1;
        x = diff(barpositioncell{n}(1:thc(n)));% based on continuous movements leftward having (-) diff
        if any(x >= 0);
            y = max(find(x >= 0)) + 1;
            initiation(n) = time_trial{n}(y) - timestart(n);
        else initiation(n) = time_trial{n}(1) - timestart(n);    % if all diff (-), first movement was in direction of threshold movement
            clear x;
            clear y;
        end
    elseif choice(n) == 1;
        x = diff(barpositioncell{n}(1:thc(n)));                  % based on continuous movements rightward having (+) diff
        if any(x <= 0)
            y = max(find(x <= 0)) + 1;
            initiation(n) = time_trial{n}(y) - timestart(n);
        else initiation(n) = time_trial{n}(1) - timestart(n);    % if all diff (+), first movement was in direction of threshold movement
            clear x;
            clear y;
        end;
    end
end



%% make 'taskbase_bar' structure

%   (:) makes all fields single column vectors
taskbase_bar.trajectory=barpositioncell(:); % all trial trajectories (cell array)
taskbase_bar.trial_start= timestart(:); % start times for all trials (Nx1 vector)
taskbase_bar.luminance= brightness_trial(:); % lum of bar for all trials (Nx1 vector)
taskbase_bar.reward_side=reward_side(:);% left = -1, right = 1 (Nx1 vector)
taskbase_bar.choice=choice(:);% actual side chosen; left = -1, right = 1, no choice = 0 (Nx1 vector)
taskbase_bar.initiation=initiation(:);% time of 'choice' movement start (Nx1 vector)
taskbase_bar.expid=rawdatafile; % gives file name
taskbase_bar.displacement_correct=th; % movement size from starting position needed to be considered a choice (1-element vector)
% taskbase_bar.displacement_bound=bound; %absolute limit position that the bar can move away from target (1-element vector)

%%