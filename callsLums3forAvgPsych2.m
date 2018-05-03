function [error, positions] = callsLums3forAvgPsych(taskbase, contra, ipsi, lumTrials)

%Updated 9/11/17 for plotting trajectories - not working yet; 
%Purpose: plot the average contra choice for psychometric function

%Run this after behaveLoad6workingOppStim_NewBin_LvR to call lums3forAvg.m
%Hard code inputs: should-mv contra v. ipsi trials on line 35-36

%% Change hard-coded lumVals here:
%Pull out relavent velocities, positions, and rxTimes:
fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2';
filenamestr = 'SC2';

%% Load taskbase structure for behavioral data 
file1 = 'mSC2lumstim_2017_08_23_16_tb.mat';
file2 = 'mSC2lumstim_2017_08_27_15_tb.mat';
file3 = 'mSC2stimlums_2017_09_01_16_tb.mat';
file4 = 'mSC2stimlum_2017_09_07_16_tb.mat';
numInputFiles = 4;

% Input luminance values - easier to find uniques, but set up the psych func's sensitivity parameters
lumVal_L1 = 255;   %should-mv left trials; ipsi for SC2
lumVal_L2 = 64;
lumVal_L3 = 15;

lumVal_R1 = -255;  %should-mv right trials; contra for SC2
lumVal_R2 = -64;
lumVal_R3 = -15;

lumVals = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L3, lumVal_L2, lumVal_L1]'; 


S = load(file1);
taskbase = S.taskbase;
trialMat = taskbase.trialMat;
timeTrial = taskbase.time_trial;
clear S;


for i = 1:length(numInputFiles)
    [yStim(i), num_L_stim(i), num_R_stim(i), yUnstim(i), num_L_unstim(i), num_R_unstim(i)] = lums3forAvg(trialMat, timeTrial, lumVals)
    yStimAll(:,i) = yStim(i);
    num_L_stimAll(:,i) = num_L_stim(i);
    num_R_stimAll(:,i) = num_R_stim(i);
    yUnstimAll(:,i) = yUnstim(i);
    num_L_unstimAll(:,i) = num_L_unstim(i);
    num_R_unstimAll(:,i) = num_R_unstim(i);
end

    %% Plotting for the stim v unstim plot:    
    figure; hold on;
    x2 = [-257 -66 -17 13 62 253]; %to jitter the unstim results so they don't plot over stim
    x2 = x2';
    plotUnstim = plot(x2, yUnstimAll, 'k.', 'MarkerSize', 18);
    hold on;
%     best_fit_Unstim = glmfit(x2, [num_R_unstim (num_R_unstim + num_L_unstim)], 'binomial');
%     x_axis_fit = -255:255;
%     y_fit = glmval(best_fit_Unstim, x_axis_fit, 'logit');
%     plot(x_axis_fit, y_fit, 'k-');

    %Stim trials
    xStim = lumVals; 
    plotStim = plot(xStim, yStimAll, 'bl.', 'MarkerSize', 18);
%     yStim = [fracCorr_R1Stim, fracCorr_R2Stim, fracCorr_R3Stim, fracIncorr_L3Stim, fracIncorr_L2Stim, fracIncorrHighLums_L1Stim]';   
%     plotStim = plot(xStim, yStim, 'bl.', 'MarkerSize', 18);
    set(gca, 'XLim', [-260 260], 'Ylim', [0 1]);
        xlabel('Luminance vals', 'FontSize', 14);
        ylabel('Fraction contra-mv choice', 'FontSize', 14);
%     num_L_stim = [numelincorrHighRmvStim, numelincorrLowRmvStim, numelincorrLowestRmvStim, numelcorrLowestLmvStim, numelcorrLowLmvStim, numelcorrHighLmvStim]';
%     num_R_stim = [numelcorrHighRmvStim, numelcorrLowRmvStim, numelcorrLowestRmvStim, numelincorrLowestLmvStim, numelincorrLowLmvStim, numelincorrHighLmvStim]';
% 
%     best_fit_stim = glmfit(xStim, [num_R_stim (num_R_stim + num_L_stim)], 'binomial');                          
%     x_axis_fit = -255:255;
%     y_fitStim = glmval(best_fit_stim, x_axis_fit, 'logit');
%     plot(x_axis_fit, y_fitStim, 'bl-');

%     mouseInfo = strcat(mouseName, ',', num2str(numSession));
%     title(mouseInfo);

return
