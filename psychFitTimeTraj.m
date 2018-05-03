
function [alignTimes, positions] = psychFitTimeTraj(fileDirectory, lumVal_L1, lumVal_L2, lumVal_L3,...
    lumVal_L4, lumVal_R1, lumVal_R2, lumVal_R3, lumVal_R4, tHold, maxNumTrials)

%Wrapper function
%Updated 5/15/15 for plotting recently refined trajectories (Vgatone and Vgattwo)
%Purpose of this script is to plot positions of bar over trial duration(time)

%This script calls PsychRxnTime3LumsPos.m thresholded positions for
%correct v. incorrect trials 11/30/14

%BS Tweaked for renewal submission 3/4/15

%% Change hard-coded lumVals here:
%//////////////////////////////////////////////////////////////////////////
%Pull out relavent velocities, positions, and rxTimes:
% cd '/Users/stubblefielde/Desktop/VGATone/easiest/150512/'; 
% fileDirectory = '/Users/stubblefielde/Desktop/VGATone/easiest/150512/'         %for all trials (4/15/15) 
fileDirectory = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks/';
cd(fileDirectory);
tHold = 50;  %this is still a velocity t-hold
maxNumTrials = 150;

%L trials   
lumVal_L1 = -255;
lumVal_L2 = -255;
lumVal_L3 = -255;
lumVal_L4 = -255;

%R trials
lumVal_R1 = 255;
lumVal_R2 = 255;
lumVal_R3 = 255;
lumVal_R4 = 255;

%cd '/Users/stubblefielde/Desktop/mice/mice2/';

% trialsToPlot1 = 10;
% trialsToPlot2 = 20; 
% trialsToPlot3 = 30; 
% 
%//////////////////////////////////////////////////////////////////////////
%fileDirectory = '/Users/stubblefielde/Desktop/VGATone/easiest/moreTrials/'    %for all trials (4/15/15) 
% cd '/Users/stubblefielde/Desktop/VGATone/easiest/150512/';                     %updated 5/15/15
% cd (fileDirectory);
maxNumTrials = 250;
% tHold = 13;    %update this 2.28.16... was the hard-coded position that is reached to denote incorrect trials

% fileDirectory = '/Users/stubblefielde/Desktop/mice/mice2/';                 %for visitor re-approval trajectories

    [trialTimeVec, trialVec, numTrials, maxTrials, mousenumber, session] = psychFitRxnTime3LumsPos3(fileDirectory);

%% WORK ON THIS
    trialNum = 1:(numTrials);    
    newTimeMat = NaN(maxTrials, 1:trialNum);                                                
    for i = 1:length(trialNum)
        newTimeMat(1:maxTrials, i) = trialTimeVec(1:maxTrials, i) - trialTimeVec(1, i);
    end                                                                         %These are new times aligned to trialStartTimes
 
    
% Sanity check for trialVec = exact initial bar position :
    initPos = trialVec(1,:)
%     for p = 1:length(initPos)
%         nineties = find(initPos(p)) == 90
%     end
    oneO5s = find(trialVec(1,:)) == abs(105)
    oneTwenties = find(trialVec(1,:)) == abs(120)
    oneThirty5s = find(trialVec(1,:)) == abs(135)
%     end

%% Now plot position of bar with respect to time
% For x-axis = time relative to the trial start time
%     figure;
%     hold on;
    plot(newTimeMat(:,trialsToPlot1), trialVec(:,trialsToPlot1), 'k', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot2), trialVec(:,trialsToPlot2), 'bl', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot3), trialVec(:,trialsToPlot3), 'r', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot4), trialVec(:,trialsToPlot4), 'm', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot5), trialVec(:,trialsToPlot5), 'g', 'LineWidth',2);
    plot(newTimeMat(:,trialsToPlot6), trialVec(:,trialsToPlot6), 'bl', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot7), trialVec(:,trialsToPlot7), 'r', 'LineWidth',2);


    
%      plot(trialVec(:,trialsToPlot1), newTimeMat(:,trialsToPlot1), 'k', 'LineWidth',2);
%      plot(trialVec(:,trialsToPlot2), newTimeMat(:,trialsToPlot2), 'g', 'LineWidth',2);
% 
%      plot(trialVec(:,trialsToPlot3), newTimeMat(:,trialsToPlot3), 'bl', 'LineWidth',2); %keep 
%      plot(trialVec(:,trialsToPlot4), newTimeMat(:,trialsToPlot4), 'k', 'LineWidth',2);
% 
%      plot(trialVec(:,trialsToPlot5), newTimeMat(:,trialsToPlot6),'r', 'LineWidth',2);
%      plot(trialVec(:,trialsToPlot6), newTimeMat(:,trialsToPlot5), 'm', 'LineWidth',2);

%     title(strcat(mousenumber,  session, '_ Unthresholded positions'));
%     legend('trial 33','trial 50', 'trial 59')   %, 'trial 40', 'trial 50', 'location', 'NorthEast');
%     set(gca, 'XLim', [0 1500]);
%     set(gca, 'XTickLabel', [0:250:1500]);
%     ylabel('Time from trial start (ms)', 'FontSize', 24);
%     xlabel('(L)        Bar position         (R)', 'FontSize', 24);
%     set(gca, 'XLim', [-300 300], 'FontSize', 20);
%     set(gca, 'Ylim', [-150 150], 'FontSize', 20);
%     hold off
% 
 %%
    clear data;
    clear LRtrialMove;
    clear time;
    clear timeDiff;
    clear timeDiffmsPerIndex;
    clear v;
    clear vals;
    clear rxnTime;
    clear data;
    clear Ltot100_lum;                        
    clear Ltot66_lum;
    clear Ltot33_lum; 
    clear LtotR_lum; 
    clear Rtot100_lum; 
    clear Rtot66_lum; 
    clear Rtot33_lum; 
    clear RtotL_lum; 
    clear LchoiceCor;
    clear RchoiceCor;
    clear inds;
    clear vel;
    clear velR;
    clear velTrace;
    clear luminance;
    clear btTrials
    clear x;
    clear y;
    clear trialNum;
    clear LcorLum;
    clear RcorLum;
    clear trialnumLcor;
    clear trialnumRcor;
    clear trialStartInds;
    clear trialEndInds;
    clear trialStartTimes;
    clear trialEndTimes;
    clear trialVec;
    clear trialt;
    clear trialx;
    clear tTrial;
    clear trialT;
    clear trialI;
    clear trialTimeVec;
    clear initCorInds;
    clear maxTrials;
    clear newCorInds;
    clear trialVec;
    clear trialsToPlot3 = 30;
    clear trialsToPlot4 = 40;
    clear trialsToPlot5 = 50;

return

