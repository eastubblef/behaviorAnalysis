function [error, positions] = callsLums3forAvgPsych(taskbase, contra, ipsi, lumTrials)

%Updated 9/14/17 for plotting trajectories - WORKS!
%Purpose: plot the average contra choice for psychometric function

%Run this after behaveLoad6workingOppStim_NewBin_LvR 
%This function calls lums3forAvg
%Hard code inputs: should-mv contra v. ipsi trials on line 35-36

%% Change hard-coded lumVals here:
%Pull out relavent velocities, positions, and rxTimes:
fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2';
filenamestr = 'SC2';
stimEpochEnd = 750;   %Default is to set this high to include all trials %for pulling out trials whose RX times occurred within the optostim epoch

%% Load taskbase structure for behavioral data 
% fileName = load('mSC2stimlum_2017_09_07_16_tb.mat');
% numInputFiles = 4;
S = load('mSC2stimlums_2017_09_01_16_tb.mat');
taskbase = S.taskbase;
trialMat = taskbase.trialMat;
trialMat12 = trialMat(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat = trialMat(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat = newTrialMat;
clear S;

S1 = load('mSC2lumstim_2017_08_23_16_tb.mat');
taskbase1 = S1.taskbase;
trialMat1 = taskbase1.trialMat;
trialMat12 = trialMat1(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat1 = trialMat1(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat1 = newTrialMat1;
clear S1;

S2 = load('mSC2lumstim_2017_08_27_15_tb.mat');
taskbase2 = S2.taskbase;
trialMat2 = taskbase1.trialMat;
trialMat12 = trialMat2(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat2 = trialMat2(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat2 = newTrialMat2;
clear S2;

S3 = load('mSC2stimlums_2_2017_08_30_15_tb.mat');
taskbase3 = S3.taskbase;
trialMat3 = taskbase3.trialMat;
trialMat12 = trialMat3(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat3 = trialMat3(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat3 = newTrialMat3;
clear S3;

S4 = load('mSC2stimlum_2017_09_07_16_tb.mat');
taskbase4 = S4.taskbase;
trialMat4 = taskbase4.trialMat;
trialMat12 = trialMat(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat4 = trialMat4(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat4 = newTrialMat4;
clear S4;

S5 = load('mSC2lumstim_2017_09_12_14_tb.mat');
taskbase5 = S5.taskbase;
trialMat5 = taskbase5.trialMat;
trialMat12 = trialMat(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat5 = trialMat5(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat5 = newTrialMat5;
clear S5;

S6 = load('mSC2lums_stim_2017_08_18_17_tb.mat');
taskbase6 = S6.taskbase;
trialMat6 = taskbase6.trialMat;
trialMat12 = trialMat(:,12);  %rx times
rxEpochTrials = find(trialMat12 <= stimEpochEnd);
newTrialMat6 = trialMat6(rxEpochTrials,:);  %new trial matrix, based on rx times within the optostim epoch
trialMat6 = newTrialMat6;
clear S6;

% Input luminance values - easier to find uniques, but set up the psych func's sensitivity parameters
lumVal_L1 = 255;   %should-mv left trials; ipsi for SC2
lumVal_L2 = 64;
lumVal_L3 = 15;

lumVal_R1 = -255;  %should-mv right trials; contra for SC2
lumVal_R2 = -64;
lumVal_R3 = -15;

lumVals = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L3, lumVal_L2, lumVal_L1]'; 

[yStim, num_L_stim, num_R_stim, yUnstim, num_L_unstim, num_R_unstim] = lums3forAvg(trialMat, lumVals)
[yStim1, num_L_stim1, num_R_stim1, yUnstim1, num_L_unstim1, num_R_unstim1] = lums3forAvg(trialMat1, lumVals)
[yStim2, num_L_stim2, num_R_stim2, yUnstim2, num_L_unstim2, num_R_unstim2] = lums3forAvg(trialMat2, lumVals)
[yStim3, num_L_stim3, num_R_stim3, yUnstim3, num_L_unstim3, num_R_unstim3] = lums3forAvg(trialMat3, lumVals)
[yStim4, num_L_stim4, num_R_stim4, yUnstim4, num_L_unstim4, num_R_unstim4] = lums3forAvg(trialMat4, lumVals)
[yStim5, num_L_stim5, num_R_stim5, yUnstim5, num_L_unstim5, num_R_unstim5] = lums3forAvg(trialMat5, lumVals)
[yStim6, num_L_stim6, num_R_stim6, yUnstim6, num_L_unstim6, num_R_unstim6] = lums3forAvg(trialMat6, lumVals)

yStimAll = horzcat(yStim, yStim1, yStim2, yStim3, yStim4, yStim5);

num_L_stimAll = horzcat(num_L_stim, num_L_stim1, num_L_stim2, num_L_stim3, num_L_stim4, num_L_stim5, num_L_stim6);
num_R_stimAll = horzcat(num_R_stim, num_R_stim1, num_R_stim2, num_L_stim3, num_R_stim4, num_R_stim5, num_R_stim6);
yUnstimAll = horzcat(yUnstim, yUnstim1, yUnstim2, yUnstim3, yUnstim4, yUnstim5, yUnstim6);
num_L_unstimAll = horzcat(num_L_unstim, num_L_unstim1, num_L_unstim2, num_L_unstim3, num_L_unstim4, num_L_unstim5, num_L_unstim6);
num_R_unstimAll = horzcat(num_R_unstim, num_R_unstim1, num_R_unstim2, num_R_unstim3, num_R_unstim4, num_R_unstim5, num_R_unstim6);

%% Plotting for the stim v unstim plot:  
% Plot ea. individual point per session
% Unstim trials:
figure; hold on;
x2 = [-257 -66 -17 13 62 253]; %to jitter the unstim results so they don't plot over stim
x2 = x2';
plotUnstim = plot(x2, yUnstimAll, 'k.', 'MarkerSize', 18);

%Stim trials
xStim = lumVals; 
plotStim = plot(xStim, yStimAll, 'bl.', 'MarkerSize', 18);
    set(gca, 'XLim', [-260 260], 'Ylim', [0 1]);
    xlabel('Luminance vals', 'FontSize', 14);
    ylabel('Fraction contra-mv choice', 'FontSize', 14);
    hold off;

%% Plotting for the stim v unstim plot:  
% average the stim & unstim points across sessions
% unstim points:
figure; hold on;

x = lumVals;
x2 = [-257 -66 -17 13 62 253]; %to jitter the unstim results so they don't plot over stim
x2 = x2';
yUnstimAll = horzcat(yUnstim, yUnstim1, yUnstim2, yUnstim3, yUnstim4, yUnstim5);
yUnstimAll_avg = mean(yUnstimAll,2);
yUnstimAll_stdev = std(yUnstimAll, 0, 2);
eUnstim = errorbar(x2, yUnstimAll_avg, yUnstimAll_stdev,'k.', 'MarkerSize', 25);
eLineUnstim = eUnstim.LineWidth;
eUnstim.LineWidth = 1.5;
% For the fit:
num_L_unstimAll = horzcat(num_L_unstim, num_L_unstim1, num_L_unstim2, num_L_unstim3, num_L_unstim4, num_L_unstim5, num_L_unstim6);
num_L_unstimAll_avg = mean(num_L_unstimAll,2);
num_R_unstimAll = horzcat(num_R_unstim, num_R_unstim1, num_R_unstim2, num_R_unstim3, num_R_unstim4, num_R_unstim5, num_R_unstim6);
num_R_unstimAll_avg = mean(num_R_unstimAll,2);
    best_fit_Unstim = glmfit(x2, [num_R_unstimAll_avg (num_R_unstimAll_avg + num_L_unstimAll_avg)], 'binomial');
    x_axis_fit = -255:255;
    y_fit = glmval(best_fit_Unstim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fit, 'k-', 'lineWidth', 2);

% stim points:  
yStimAll_avg = mean(yStimAll, 2);
yStimAll_stdev = std(yStimAll, 0, 2);
eStim = errorbar(x, yStimAll_avg, yStimAll_stdev, 'bl.', 'MarkerSize', 25);
eLineStim = eStim.LineWidth;
eStim.LineWidth = 1.5;
%For the fit:
num_L_stimAll = horzcat(num_L_stim, num_L_stim1, num_L_stim2, num_L_stim3, num_L_stim4, num_L_stim5, num_L_stim6);
num_L_stimAll_avg = mean(num_L_stimAll,2);
num_R_stimAll = horzcat(num_R_stim, num_R_stim1, num_R_stim2, num_L_stim3, num_R_stim4, num_R_stim5, num_R_stim6);
num_R_stimAll_avg = mean(num_R_stimAll,2);
    best_fit_stim = glmfit(xStim, [num_R_stimAll_avg (num_R_stimAll_avg + num_L_stimAll_avg)], 'binomial');                          
    x_axis_fit = -255:255;
    y_fitStim = glmval(best_fit_stim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fitStim, 'bl-', 'lineWidth', 2);

    set(gca, 'XLim', [-260 260], 'Ylim', [0 1]);
    xlabel('Luminance vals', 'FontSize', 16);
    ylabel('Fraction contra-mv choice', 'FontSize', 16);
    hold on;

%     mouseInfo = strcat(mouseName, ',', num2str(numSession));
     title('6 stim sessions');

return
