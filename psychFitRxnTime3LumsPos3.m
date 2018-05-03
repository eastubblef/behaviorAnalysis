
function [trialTimeVec, trialVec, numTrials, maxTrials, mousenumber, session] = psychFitRxnTime3LumsPos3(fileDirectory, tHold, maxNumTrials)

%USE W/ psychFitTimeTraj.m WRAPPER

%Created BS (DudmanLab) 10/5/14: 
%Updated BS 10/12/14 to provide more meaningful Rx times (t-hold on velocity)
%Updated BS 10/13/14 to provide movement trajectories over trial duration (NOTE: no t-hold on velocity)
%Updated BS 12/10/14 so Gidon could make sense of this
%Updated BS 2/20/15
%Updated BS 4/20/15 for full 100% lums and closer Gaussian to center
%Updated BS 2.28.16 for blinking lum task & 2AFC

% Purpose of this function is to generate rx time plots for more difficult luminance task
% Additionally, Gidon wants raw trajectories (position of bar, for now, across ea. trial)

% Inputs: 1 file/session at a time is called within a loop and plotted from
%         wrapper called psychFitTime3_pos.m

% Extracted data structure 1 = data.trajectories is more "raw" 
%      time = data.trajectories(:, 1);
%      x = data.trajectories(:, 2);
%      y = data.trajectories(:, 3);
%      btTrials = data.trajectories(:, 5);
%      luminance = data.trajectories(:, 6);  

% Extracted data structure 2 = data.trials is more refined
%     data.trials(:,1) = trialNum;
%     data.trials(:,2) = trialStartInds;
%     data.trials(:,3) = trialEndInds;
%     data.trials(:,4) = trialStartTimes;
%     data.trials(:,5) = trialEndTimes;
%     data.trials(:,6) = LRtrial;
%     data.trials(:,7) = LchoiceCor;
%     data.trials(:,8) = RchoiceCor;
%     data.trials(:,9) = lumTrials;
    
% Inputs: 
%   fileDirectory for generating extracted information for all behavioral sessions within that folder
%   lumVals are the decided-upon luminance vals to make the task more difficult
%   tHold indicates the value of velocity change (from zero) to pull out meaningful Rx times
%   maxNumTrials needed here since Java state at maxNumTrials+1 during aquisition renders no pick-up 

% Outputs: 
%   structure called "data" in which data.trials field contains relavent information
%   plot of raw positions (trajectories) over trial duration per session
%   subplot of t-holded velocity-related rx times for same session

%% These vals should be worked into the wrapper - for now, they appear below in ll. 310-316
%////////////////////
% trialsToPlot1 = 10;
% trialsToPlot2 = 20;
% trialsToPlot3 = 30;
%///////////////////

% fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks'
% filenamestr = 'mVgatthreetag';
fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16'
filenamestr = 'SC2';

    if exist('filenamestr', 'var');
        [filenamestrE, path] = uigetfile('*.csv*','select the .csv file', fileDirectory)               %BS get the actual # correct trials (and TS) from csv file
    end
    
%     if exist('filenamestr', 'var');
%         [filenamestrF, path] = uigetfile('*_p.csv*','select the _p.csv file', fileDirectory)            %pxy file has solenoid discharge TS with respect to trial starts
%     end
    
    if exist('filenamestr', 'var');
        [filenamestrX, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fileDirectory)            %pxy file has solenoid discharge TS with respect to trial starts
    end
    
    %% Load csv behavior data
    
    clear disp* trial* 
     skipInitTrials = 0;

     csvFile = filenamestrE;
     data.csvTrials = dlmread(csvFile,',',2,0);
     data.attemptedTrials = data.csvTrials(:,1);
     data.numRewardedTrials = data.csvTrials(end,1);
     data.trialStarts = data.csvTrials(:,2);
     data.trialEnds = data.csvTrials(:,4);

     pxyFile = filenamestrX;
     data.pxyTS = dlmread(pxyFile,',',2,0);  % can use the 'taskbase' file here to parse out solenoid TS & centered TS
     data.trajectories = data.pxyTS;
%      pFile = filenamestrF;
%      data.trialParams = dlmread(pFile,',',1+skipInitTrials,0);        %11th col is the column containing information about which trials were stimulated     
    
     time = data.trajectories(:, 1);
     x = data.trajectories(:, 2);
     y = data.trajectories(:, 3);
     btTrials = data.trajectories(:, 5);
     luminance = data.trajectories(:, 6);    %for blinking
    
    %% General information about the file //BS Pull out relevant info 8/29/14 - DON'T USE YET (12/10/14)
%      % Note that these vals within this block will have fewer indices than vars to follow
%      
%      indArray = 2:size(data.trialTimes,1);                                                  % From the .csv file
%      data.startToEvent  = data.trialTimes(indArray,3) - data.trialTimes(indArray,2);        % use? time from trial start to threshold cross (lever)
%      data.endToStart    = data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4);      % inter-trial interval (should be 3 sec)
%      data.iti = (data.trialTimes(indArray,2) - data.trialTimes(indArray-1,4));              % From the .csv file...I'd rather define it as ITI 
%     
%      data.eventToEnd    = data.trialTimes(indArray,4) - data.trialTimes(indArray,3);        % use? reward delay
%      numTrials   = size(data.trialTimes,1);                                                 % number of total trials
%      lastTrial = data.trialTimes(end,1);                                                    % will be offset by 2 (defined above)
%  
%      %BS IMPORTANT: Shift var represents the alignment of continuous data with the trial times data (eg. state of 0 becomes 2, in the last col, and that first col num becomes the trial-start time
%      %shift = data.trajectories.contXY(find(data.trajectories.contXY(:,4)>0, 1,'first'),1) - data.trialTimes(1,2); % calculate shift

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
    vel = abs(x);                                                           %pick up any movement
%    velThold = find(vel > tHold);  NOT IN USE                              %threshold velocity, lower bound
    
    trialStartInds = find(diff(btTrials) == -1) +1;                         %transition from 1 to 0; trial start 
    trialStartTimes = time(trialStartInds); 
    
    %Updated 2.28.16 to find proper initial positions
    %Find the initial bar positions from pxy file: will be one interative value before the TS for trialStart
    forPosTrialStarts = trialStartInds-1;
    posTrialStarts = y(forPosTrialStarts);
    
    trialEndInds = find(diff(btTrials) == 1) +1;                            %these are the last 0s before 1 (before next ITI)
    trialEndTimes = time(trialEndInds);
    
    if numel(trialStartInds) > numel(trialEndInds)
        trialStartInds = trialStartInds(1:end-1);                           %adapt for starting a trial that never ended
    end

    if numel(trialStartTimes) > numel(trialEndTimes)
        trialStartTimes = trialStartTimes(1:end-1);                         %adapt for starting a trial that never ended
    end   
    
    if numel(trialStartTimes) == maxNumTrials + 1 && numel(trialStartInds) == maxNumTrials + 1
        trialStartTimes = trialStartTimes(1:end-1);
        trialStartInds = trialStartInds(1:end-1);
    end
    
    if numel(trialEndTimes) == maxNumTrials + 1 && numel(trialEndInds) == maxNumTrials + 1
        trialEndTimes = trialEndTimes(1:end-1);
        trialEndInds = trialEndInds(1:end-1);
    end
    
    %% Assess trial R v L and correct choices
    clear LchoiceCor;
    clear RchoiceCor;
    clear LRtrialMove;
    clear newCorInds;
    
    velR = x;
    numSamps = numel(vel);
    numSampsToUse = 300;
%     velTrace = zeros(numel(trialStartInds), numSampsToUse);
    
    trialNum = 1:length(trialStartInds);
      lumR = unique(luminance);                                             %these just represent the individual lums used
      lumL = -unique(luminance);
      lumAll = vertcat(lumR, lumL);
      lumTrials = luminance(trialStartInds);                                %vector of luminance vals per trial
          
    LRtrial = y(trialStartInds);                                            %these are the initial bar/wheel positions
    
    initCorInds = zeros(size(trialStartInds)); 
    for t = 1:length(trialStartInds)
       % initCorInds(t) = find(vel(trialStartInds(t):numSamps) >1 ,1)+ 1;    %will be a vector of the index# of velocity change for any movement, with respect to trialStartInds
        initCorInds(t) = find(vel(trialStartInds(t):numSamps) > tHold,1);    %will be a vector of the trialStartInd+initCorInd for movements with vel > tHold
       % initCorInds(t) = y(max(vel(trialStartInds(t):trialEndInds(t),1)));
    end
    initCorInds = initCorInds';                                             
    LRtrialMove = y(initCorInds);                                           %indexing this way is not accurate since trialStartTimes/Inds are not factored in
    
    %Turn this into an interpretable vector - rxn times also calculated from this (l. 530)
    for new = 1:length(initCorInds)
        newCorInds(new) = initCorInds(new) + trialStartInds(new);
    end
    newCorInds = newCorInds';
    
    LRtrialMove = y(newCorInds);                                            %value of that movement @ t-holded velocity
     
    %L and R correct choices
    for l = 1:length(LRtrial)
        Ltrial = find(LRtrial < 1);
        LchoiceCor(l) = LRtrial(l) < 1 & LRtrialMove(l) > LRtrial(l);
        
        Rtrial = find(LRtrial > 1);
        RchoiceCor(l) = LRtrial(l) > 1 & LRtrialMove(l) < LRtrial(l);
    end
    LchoiceCor = LchoiceCor';
    RchoiceCor = RchoiceCor';
    
    %% New matrix of extracted info = data.trials field 
    if exist('data');
       clear data;
    end
    data.trials(:,1) = trialNum;
    data.trials(:,2) = trialStartInds;
    data.trials(:,3) = trialEndInds;
    data.trials(:,4) = trialStartTimes;
    data.trials(:,5) = trialEndTimes;
    data.trials(:,6) = LRtrial;
    data.trials(:,7) = LchoiceCor;
    data.trials(:,8) = RchoiceCor;
    
    %Get luminances for L v R trials 
    for p = 1:length(LRtrial)
        if LRtrial(p) < 1
            lumTrials(p) = -1 * lumTrials(p);
        end
    end
    
    data.trials(:,9) = lumTrials;
    
    %% Find positions (and times) from trialStart to trialEnd for Gidon

    maxTrials = max(trialEndInds - trialStartInds);    
    trialVec = NaN(maxTrials, length(trialNum));
    trialTimeVec = NaN(maxTrials, length(trialNum));                        %To fit in all times in trialTimeVec
    trialI = abs(trialEndInds - trialStartInds);
    
    %Use time to index position:
    timeDiff = trialEndTimes - trialStartTimes;
    timeDiffmsPerIndex = timeDiff./trialI;
    timesNew = cumsum(timeDiffmsPerIndex);
    
    data.trials(:,10) = timeDiffmsPerIndex;                                 %In case I need the average time per index of that trial 
    data.trials(:,11) = timesNew;                                           %In case I need the average cummulative sum of the ms/trial index
    
    trialT = timeDiffmsPerIndex .* (trialI);    
    
    for iTrial = 1:length(trialI)
        trialVec(1:trialI(iTrial), iTrial) = y(trialStartInds(iTrial):trialEndInds(iTrial)-1);          %These are the trial indices with position
    
        trialTimeVec(1:trialI(iTrial), iTrial) = time(trialStartInds(iTrial):trialEndInds(iTrial)-1);   %These are raw times (they don't start at zero) to plot with respect to position
    end
    
    %% Plot bar position across ea. trial (indices) - 1st figure
     % x-axis = trial indices from trialStart to trialEnd
    
    figure;
    hold on;
    plot(trialVec(:,1), 'k', 'LineWidth',2);
    plot(trialVec(:,3), 'bl', 'LineWidth',2);
%     plot(trialVec(:,40), 'm', 'LineWidth',2);
%     plot(trialVec(:,50), 'g', 'LineWidth',2);
    
    title(strcat(mousenumber,  session, '_ Unthresholded positions'));
    legend('trial 10','trial 20', 'trial 30', 'trial 40', 'trial 50', 'location', 'NorthEast');
    
%     if maxTrials > 1500
%         set(gca, 'XLim', [0 1500]);
%         else set(gca, 'XLim', [0 (maxTrials/2)]);
%     end
    
    xlabel('Trial Indices', 'FontSize', 16);
    ylabel('Position', 'FontSize', 16);
    set(gca, 'YLim', [-200 200], 'FontSize', 15);
    set(gca, 'XLim', [0 500], 'FontSize', 15);
     hold off;
        
%%  Now plot position of bar with respect to time - 2nd figure
    %For x-axis = time relative to the trial start time
    %CLEAN THIS UP TO RETURN TO WRAPPER EVENTUALLY
    
    trialNum = 1:(numTrials);    
    newTimeMat = NaN(maxTrials, 1:trialNum);                                                
    for i = 1:length(trialNum)
        newTimeMat(1:maxTrials, i) = trialTimeVec(1:maxTrials, i) - trialTimeVec(1, i); %These are new times aligned to trialStartTimes
    end
    
%     trialsToPlot1 = 10;
%     trialsToPlot2 = 20;
%     trialsToPlot3 = 30;
%     trialsToPlot4 = 40;
%     trialsToPlot5 = 50;
%     trialsToPlot6 = 60;
%     trialsToPlot7 = 70; 

    trialsToPlot1 = 1;
    trialsToPlot2 = 2;
    trialsToPlot3 = 3;
    trialsToPlot4 = 4;
    trialsToPlot5 = 5;
    trialsToPlot6 = 6;

    figure;
    hold on;
    plot(newTimeMat(:,trialsToPlot1), trialVec(:,trialsToPlot1), 'k', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot2), trialVec(:,trialsToPlot2), 'bl', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot3), trialVec(:,trialsToPlot3), 'r', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot4), trialVec(:,trialsToPlot4), 'm', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot5), trialVec(:,trialsToPlot5), 'g', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot6), trialVec(:,trialsToPlot6), 'k', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot7), trialVec(:,trialsToPlot7), 'bl', 'LineWidth',2);

    
    xlabel('Time from trial start (ms)', 'FontSize', 16);
    ylabel('Position', 'FontSize', 16);
%      set(gca, 'XLim', [0 10000], 'FontSize', 15);
    set(gca, 'XLim', [0 6000], 'FontSize', 15);
    set(gca, 'YLim', [-350 350], 'FontSize', 15);
    title(strcat(mousenumber,  session, '_ Unthresholded positions'));
    legend('trial 10','trial 20', 'trial 30', 'trial 40', 'trial 50', 'trial 60', 'trial 70');
    hold off;
    
    %% Determine % correct for each luminance - 12/10/14 USE RATHER PSYCH4LUMSLOW.M FOR PSYCH F(X) PLOTS
    
    Ltot100_lum = find(lumTrials == lumVal_L1);                             %update these lums for sessions of diff # lum vals; see first lines
    Ltot66_lum = find(lumTrials == lumVal_L2);
    Ltot33_lum = find(lumTrials == lumVal_L3);
    LtotR_lum = find(lumTrials == lumVal_L4);
    
    Rtot100_lum = find(lumTrials == lumVal_R1);
    Rtot66_lum = find(lumTrials == lumVal_R2);
    Rtot33_lum = find(lumTrials == lumVal_R3);
    RtotL_lum = find(lumTrials == lumVal_R4);
        
    %For the adjusted x axis (so no gaps remain)
%     lumVal_L1 = -4; lumVal_L2 = -3; lumVal_L3 = -2; lumVal_L4 = -1; lumVal_R4 = 1, lumVal_R3 = 2, lumVal_R2 = 3, lumVal_R1 = 4;
%     lumRightAdjust = ['-4', '-3', '-2', '-1', '0', '1', '2', '3'];

    for al = 1:length(lumTrials)
        if lumTrials(al) == lumVal_L1
            lumAdjust(al) = -4;
        else
            if lumTrials(al) == lumVal_L2
                lumAdjust(al) = -3;
            else
                if lumTrials(al) == lumVal_L3
                    lumAdjust(al) = -2;
                else
                    if lumTrials(al) == lumVal_L4
                        lumAdjust(al) = -1;
                    end
                end
            end
        end
    end
                  
    for ar = 1:length(lumTrials)
        if lumTrials(ar) == lumVal_R1
            lumAdjust(ar) = 3;
        else
            if lumTrials(ar) == lumVal_R2
                lumAdjust(ar) = 2;
            else
                if lumTrials(ar) == lumVal_R3
                    lumAdjust(ar) = 1;
                else
                    if lumTrials(ar) == lumVal_R4
                        lumAdjust(ar) = 0;
                    end
                end
            end
        end
    end            
    lumAdjust = lumAdjust';
    
    %Left trials correct
    for a = 1:length(Ltot100_lum)                                           %for percent correct for 100% (brightest) contrast for L trials
        corLtot100(a) = LchoiceCor(Ltot100_lum(a));
    end
    percL100 = find(corLtot100 == 1); percLCor100 = 100 * (numel(percL100)/numel(corLtot100));
    percLCor100 = 100 - percLCor100;                                        %this generates "%R choice" Y-axis
    fracLBright = numel(percL100)/numel(corLtot100);                        %fraction of choices to the correct side (in case I don't want to plot %)
    fracBright = 1 - fracLBright;

    for b = 1:length(Ltot66_lum)                                            %for percent correct for 66% contrast for L trials
        corLtot66(b) = LchoiceCor(Ltot66_lum(b));
    end
    percL66 = find(corLtot66 == 1); percLCor66 = 100 * (numel(percL66)/numel(corLtot66));
    percLCor66 = 100 - percLCor66;                                          %this generates "%R choice" Y-axis
    fracL2Bright = numel(percL66)/numel(corLtot66)
    frac2Bright = 1 - fracL2Bright;

    for c = 1:length(Ltot33_lum)                                            %for percent correct for 33% contrast for L trials
        corLtot33(c) = LchoiceCor(Ltot33_lum(c));
    end
    percL33 = find(corLtot33 == 1); percLCor33 = 100 * (numel(percL33)/numel(corLtot33));
    percLCor33 = 100 - percLCor33;
    fracL2Hard = numel(percL33)/numel(corLtot33);
    frac2Hard = 1 - fracL2Hard;

    for d = 1:length(LtotR_lum)                                             %for percent correct for no contrast for L trials
        corLtotR(d) = LchoiceCor(LtotR_lum(d));
    end
    percLR = find(corLtotR == 1); percLCorR = 100 * (numel(percLR)/numel(corLtotR));
    percLCorR = 100 - percLCorR;
    fracLHard = numel(percLR)/numel(corLtotR)
    fracHard = 1 - fracLHard;

    %Right trials correct
    for q = 1:length(Rtot100_lum)                                           %for percent correct for 100% (brightest) contrast for R trials
        corRtot100(q) = RchoiceCor(Rtot100_lum(q));
    end
    percR100 = find(corRtot100 == 1); percRCor100 = 100 * (numel(percR100)/numel(corRtot100));
    fracRBright = numel(percR100)/numel(corRtot100);

    for r = 1:length(Rtot66_lum)                                            %for percent correct for 66% contrast for R trials
        corRtot66(r) = RchoiceCor(Rtot66_lum(r));
    end
    percR66 = find(corRtot66 == 1); percRCor66 = 100 * (numel(percR66)/numel(corRtot66));
    fracR2Bright = numel(percR66)/numel(corRtot66);

    for s = 1:length(Rtot33_lum)                                            %for percent correct for 33% contrast for R trials
        corRtot33(s) = RchoiceCor(Rtot33_lum(s));
    end
    percR33 = find(corRtot33 == 1); percRCor33 = 100 * (numel(percR33)/numel(corRtot33));
    fracR2Hard = numel(percR33)/numel(corRtot33);

    for t = 1:length(RtotL_lum)                                             %for percent correct for (no) contrast for R trials
        corRtotL(t) = RchoiceCor(RtotL_lum(t));
    end
    percRL = find(corRtotL == 1); percRLCor = 100 * (numel(percRL)/numel(corRtotL));
    fracRHard = numel(percRL)/numel(corRtotL);
   
 %% Plot % R choice with psychometric fit - IGNORE FOR NOW
%     figure; 
%     hold on; 
%     lumRight = [lumVal_L1, lumVal_L2, lumVal_L3, lumVal_L4, lumVal_R4, lumVal_R3, lumVal_R2, lumVal_R1];
%     lumRightAdjust = [-4, -3, -2, -1, 0, 1, 2, 3];
% 
%     yVals = [fracBright, frac2Bright, frac2Hard, fracHard, fracRHard, fracR2Hard, fracR2Bright, fracRBright];   %to plot y axis as fraction R choice
% %     yVals = [percLCor100, percLCor66, percLCor33, percLCorR, percRLCor, percRCor33, percRCor66, percRCor100];
% 
%     z = plot(lumRightAdjust, yVals, 'k.', 'MarkerSize', 24);
%     title(strcat(month, ' ', day, ' ', mousenumber, session));
% 
%   %% For the fit:
% %   
% %     for n = 1:length(yVals)
% %         newN(n) = 1;
% %     end
%     
% %     newN = newN';
% %     lumRight = lumRight';
% %     yVals = yVals';
% %     yVals2 = [numel(corLtot100), numel(corLtot66), numel(corLtot33), numel(corLtotR), numel(corRtotL), numel(corRtot33), numel(corRtot66), numel(corRtot100)];
% % %     newY = horzcat(yVals', yVals2');
% % %     bf = glmfit(lumRight, [newY, newN], 'binomial'); %why does this not work?
% %     en = yVals2';
% %     bf = glmfit(lumRightAdjust, [yVals, en], 'binomial');
% %     
% %     x_axis = lumRightAdjust;
% %     y_fit = glmval(bf, x_axis, 'logit');
% % 
% %     % return the fit y and the x it was based on, and the calculated bias
% %     fit_values.x_axis = x_axis;
% %     fit_values.y_fit = y_fit;
% % %     fit_values.bias = 50 + (bf(1)/bf(2)); % Bobby's calculation, based on Hatim's email explanation
% % 
% %     pFit = plot(fit_values.x_axis, fit_values.y_fit, 'k');
% %     set(pFit, 'LineWidth', 2);
%    
%     set(gca, 'XTickLabel', [lumRight]);
%     xlabel('Luminance vals', 'FontSize', 16);
%     ylabel('Fraction R choice (init. movement)', 'FontSize', 16);
%     set(gca, 'YTickLabel', [0:0.1:1], 'FontSize', 15);
%     set(gca, 'YLim', [0 1]);
%     hold off;
%     
%     %% Use velocity indexed to trialStartInds for plotting Rx times - use for trial difficulty, though 10/17/14
% %     
% %     figure;
% %     hold on;
%     velR = x;
%     numSamps = numel(vel);
%     numSampsToUse = 300;
%     velTrace = zeros(numel(trialStartInds)-1, numSampsToUse);
%     
%     %This will be one trial shorter, but necessary to index properly (?)
%     for v2 = 1:length(trialStartInds)-1
%         velTrace(v2,1:numSampsToUse) = velR(trialStartInds(v2):trialStartInds(v2)+numSampsToUse-1,1)';
%     end
% %     
% % %     for pl = 1:numel(trialStartInds)-1
% % %         plotVel(pl) = plot(velTrace(pl,:));
% % %     end
% %     
% %     plot(velTrace(10,:), 'r');
% %     plot(velTrace(20,:), 'k');
% %     plot(velTrace(30,:), 'g');
% %     plot(velTrace(40,:), 'bl');
% %     plot(velTrace(50,:), 'r');
% %     plot(velTrace(60,:), 'g');
% %     
% %     title(strcat(mousenumber,  session));
% 
% %     set(gca, 'XLim', [0 numSampsToUse], 'FontSize', 15);
% %     set(gca, 'YLim', [-50 50]);
% %     ylabel('Velocity ( - left;  + right)', 'FontSize', 16); 
% %     xlabel('NumSampsToUse/trial', 'FontSize', 16);
% %     title(strcat(mousenumber,  session));
% 
% %     for v = 1:numel(trialStartInds)                                         
% % %         rxnTime(v) = find(vel(trialStartInds(v):numSamps)>1,1);             %find where velocity is not zero, at the first value
% %         rxnTime(v) = trialStartTimes(newCorInds(v));
% %         if numel(rxnTime(v)) > numel(trialStartInds(v));
% %             break
% %         else if numel(rxnTime(v)) < numel(trialStartInds(v));
% %                 break
% %             end
% %          end
% %     end
%     
%%  Rxn time plots - Figure 3 subplots
%     clear rxnTime;
%     clear rxnTimePre;
%     
%     rxnTimePre = time(newCorInds);
%     rxnTime = rxnTimePre - trialStartTimes;
%     data.trials(:,12) = rxnTime;
%     maxRx = max(rxnTime);
%     
%     %Plot all (no discrepancy bt L v R) for correct Rx times across all trials
%     figure;
%     subplot(2,2,1); hold on;
%     plot(trialNum, rxnTime, 'k.', 'MarkerSize', 16);
%     ylabel('Rxn time (ms)', 'FontSize', 14); xlabel('Trial number', 'FontSize', 14); 
%     set(gca, 'YLim', [0  maxRx]);
%     title(strcat(mousenumber,  session));
% 
%       %Plot veltrace and the correct rxn times 
% %     subplot(2,2,2);
% %     hold on;
% %    
% %     [vals,inds] = sort(rxnTime(1:end-1));                                            
% %         
% %     [mapName] = TNC_CreateRBColormap(12,'mbr');
% %     imagesc(velTrace(inds(1:end-1),:),[-100 100]); colormap(mapName);
% %     ylabel('Sorted Trials', 'FontSize', 14); xlabel('Time (samples | 19.3 ms/sample', 'FontSize', 14);
% %     plot(rxnTime(inds),1:numel(inds),'ko');
%     
%     %Plot Rx time for left v. right trials across correct trials (all lums)
%     subplot(2,2,3);
%     hold on;
%     trialnumLcor = find(LchoiceCor == 1);
%     rxnLcor = rxnTime(trialnumLcor);
%     trialnumRcor = find(RchoiceCor == 1);
%     rxnRcor = rxnTime(trialnumRcor);
%     scatL = plot(trialnumLcor, rxnLcor, 'k.', 'MarkerSize', 16);
%     scatR = plot(trialnumRcor, rxnRcor, 'r.', 'MarkerSize', 16);
%     ylabel('Rxn time (ms), correct', 'FontSize', 14); xlabel('Trial number', 'FontSize', 14); 
%     set(gca, 'YLim', [0 10000]);
%     legend('L trials correct','R trials correct', 'location', 'NorthWest');
% 
%     %Plot Rx time for luminance vals
%     subplot(2,2,4);
%     hold on;
%     LcorLum = lumAdjust(trialnumLcor);
%     RcorLum = lumAdjust(trialnumRcor); 
%     scatLlum = plot(LcorLum, rxnLcor, 'k.', 'MarkerSize', 16);
%     scatRlum = plot(RcorLum, rxnRcor, 'r.', 'MarkerSize', 16);
% %     set (gca, 'XTickLabel', [-4, -3, -2, -1, 1, 2, 3, 4]); even though these vals correspond to lum vals L(-) to R(+)
%     xlabel('Luminance vals', 'FontSize', 14);
%     ylabel('Rxn times (ms), correct', 'FontSize', 14);
%     set(gca, 'YLim', [0 10000]);
%        
% end
%  
 end