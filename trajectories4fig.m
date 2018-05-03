
function [alignTimes, positions] = trajectories(taskbase, contra, ipsi)

%Updated 8/9/17 for plotting trajectories 
%Purpose: plot positions of bar over trial duration

%Run this after behaveLoad6workingOppStim_NewBin_LvR 
%Hard code inputs: should-mv contra v. ipsi trials on line 35-36
%Hard code inputs: corrects (>,< 0) for bar centering, depending on hemisphere of optic fiber

%Examples for blks: SC1 = 9.07.17 or 9.14(?)
%Examples for rand: SC2 = 8.31 

%% Change hard-coded lumVals here:
%Pull out relavent velocities, positions, and rxTimes:
fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC1';
filenamestr = 'SC1';
% lums = 1;              %input a zero for consistent/high luminance

%% Load taskbase structure for behavioral data 
                      
if exist('filenamestr', 'var');
    [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       
end
taskbase = strcat(path, filenamestrT);
load(taskbase);

mouseName = filenamestr;
numSession = 170914;

trialMat = taskbase.trialMat;
barpositionCell = taskbase.barpositionCell;
timeTrial = taskbase.time_trial;

trialMat1 = trialMat(:,1);
trialMatStim = trialMat(:,13);
trialMatStartPos = trialMat(:,14);

%% separate out should-mv ipsi v contra trials: also change ll. 111-225 corrects
% Look at should-mv stim trials first:

% contra = find(trialMatStartPos < 1);        %for SC2; L implant - must move R to center bar
% ipsi = find(trialMatStartPos > 1);          %for SC2; L implant - must move L to center bar
% incorrContraT_hold = -400;                  %for SC2; L implant - must move R to center bar
% incorrIpsiT_hold = 400;                     %for SC2; L implant - must move R to center bar
contra = find(trialMatStartPos > 1);        %for SC1; R implant - must move L to center bar
ipsi = find(trialMatStartPos < 1);          %for SC1; R implant - must move R to center bar
incorrContraT_hold = 400;                   %for SC1; L implant - must move R to center bar
incorrIpsiT_hold = -400;                    %for SC1; L implant - must move R to center bar

contraStimAll = trialMatStim(contra);
contraStim1 = find(contraStimAll == 2);
contraStim = trialMat1(contra(contraStim1));

ipsiStimAll = trialMatStim(ipsi);
ipsiStim1 = find(ipsiStimAll == 2);
ipsiStim = trialMat1(ipsi(ipsiStim1));

% Look at should-mv no stim trials:
contraNoStim1 = find(contraStimAll == 0);
contraNoStim = trialMat(contra(contraNoStim1));
ipsiNoStim1 = find(ipsiStimAll == 0);
ipsiNoStim = trialMat1(ipsi(ipsiNoStim1));

% Turn barpositionCell into a workable matrix: ea. col is a trial's trajectories
barpositionMat = NaN(200, length(barpositionCell));
for z=1:numel(barpositionCell);
    for i=1:length(barpositionCell{z})
        barpositionMat(i,z) = barpositionCell{z}(i);
    end
end

contraNoStimTrajMat = barpositionMat(1:end-1,contraNoStim); %just get rid of the last nan
contraStimTrajMat = barpositionMat(1:end-1,contraStim);
ipsiNoStimTrajMat = barpositionMat(1:end-1,ipsiNoStim);
ipsiStimTrajMat = barpositionMat(1:end-1,ipsiStim);

%% also index into the proper times that ea. position occurred:
timeTrialMat = NaN(200, length(timeTrial));
for z=1:numel(timeTrial);
    for i=1:length(timeTrial{z})
        timeTrialMat(i,z) = timeTrial{z}(i);
    end
end

% %or, use a generic timestamp for ea. row of the barpositionMat
% timeVec = 1:length(barpositionMat);
% timeVec = 1:length(timeTrialMat(:,1))+1;
% timeVec = timeVec';

% Use the true timestamps for position withing the barpositionMat in ms
contraNoStimTrueTime = timeTrialMat(:,contraNoStim);
contraStimTrueTime = timeTrialMat(:,contraStim);
ipsiNoStimTrueTime = timeTrialMat(:,ipsiNoStim);
ipsiStimTrueTime = timeTrialMat(:,ipsiStim);

% relTime1 = timeTrialMat(:,1)-timeTrialMat(1,1);
relTime = [];
for i = 1:length(timeTrialMat(1,:))  % #cols
    relTimeMat(:,i) = timeTrialMat(:,i) - timeTrialMat(1,i);
end

contraNoStimRelTime = relTimeMat(:,contraNoStim);
contraStimRelTime = relTimeMat(:,contraStim);
ipsiNoStimRelTime = relTimeMat(:,ipsiNoStim);
ipsiStimRelTime = relTimeMat(:,ipsiStim);

%% Separate out correct should-go contras:
%find the cols for no-stim contras that had incorrects (0)
corrContraNoStim = [];
numelColsContraNoStim = numel(contraNoStimTrajMat(1,:));
corrContraNoStimCol = [];
for i = 1:length(contraNoStimTrajMat)                            %for ea row
    for j = 1:numelColsContraNoStim                              %for ea col
%         if contraNoStimTrajMat(i,j) > 0                      %SC2
        if contraNoStimTrajMat(i,j) < 0                        %SC1
            corrContraNoStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
corrContraNoStim1 = corrContraNoStimCol > 1;
corrContraNoStim = corrContraNoStimCol(corrContraNoStim1);     %these are the cols of incorr contra trials

%find the cols for stim contras that had corrects (0)
corrContraStim = [];
numelColsContraStim = numel(contraStimTrajMat(1,:));
corrContraStimCol = [];
for i = 1:length(contraStimTrajMat)                            %for ea row
    for j = 1:numelColsContraStim                              %for ea col
%         if contraStimTrajMat(i,j) > 0
        if contraStimTrajMat(i,j) < 0
            corrContraStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
corrContraStim1 = corrContraStimCol > 1;
corrContraStim = corrContraStimCol(corrContraStim1);     %these are the cols of incorr contra trials

%% Separate out incorrect should-go contras:
%find the cols for no-stim contras that had incorrects (-400 for SC2)
incorrContraNoStim = [];
numelColsContraNoStim = numel(contraNoStimTrajMat(1,:));
incorrContraNoStimCol = [];
for i = 1:length(contraNoStimTrajMat)                            %for ea row
    for j = 1:numelColsContraNoStim                              %for ea col
        if contraNoStimTrajMat(i,j) == incorrContraT_hold
            incorrContraNoStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
incorrContraNoStim1 = incorrContraNoStimCol > 1;
incorrContraNoStim = incorrContraNoStimCol(incorrContraNoStim1);     %these are the cols of incorr contra trials

%find the cols for stim contras that had incorrects (-400 for SC2)
incorrContraStim = [];
numelColsContraStim = numel(contraStimTrajMat(1,:));
incorrContraStimCol = [];
for i = 1:length(contraStimTrajMat)                            %for ea row
    for j = 1:numelColsContraStim                              %for ea col
        if contraStimTrajMat(i,j) == incorrContraT_hold
            incorrContraStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
incorrContraStim1 = incorrContraStimCol > 1;
incorrContraStim = incorrContraStimCol(incorrContraStim1);     %these are the cols of incorr contra trials

%% Separate out correct should-go ipsis:
%find the cols for no-stim contras that had corrects (0)
corrIpsiNoStim = [];
numelColsIpsiNoStim = numel(ipsiNoStimTrajMat(1,:));
corrIpsiNoStimCol = [];
for i = 1:length(ipsiNoStimTrajMat)                            %for ea row
    for j = 1:numelColsIpsiNoStim                              %for ea col
%         if ipsiNoStimTrajMat(i,j) < 0   %SC2
        if ipsiNoStimTrajMat(i,j) > 0     %SC1

            corrIpsiNoStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
corrIpsiNoStim1 = corrIpsiNoStimCol > 1;
corrIpsiNoStim = corrIpsiNoStimCol(corrIpsiNoStim1);     %these are the cols of incorr contra trials

%find the cols for stim contras that had corrects (0)
corrIpsiStim = [];
numelColsIpsiStim = numel(ipsiStimTrajMat(1,:));
corrIpsiStimCol = [];
for i = 1:length(ipsiStimTrajMat)                            %for ea row
    for j = 1:numelColsIpsiStim                              %for ea col
%         if ipsiStimTrajMat(i,j) < 0     %SC2
           if ipsiStimTrajMat(i,j) > 0    %SC1

            corrIpsiStimCol(j) = j;                        %tell me the col where corr t-hold is reached
        end
    end
end
corrIpsiStim1 = corrIpsiStimCol > 1;
corrIpsiStim = corrIpsiStimCol(corrIpsiStim1);         %these are the cols of corr contra trials

%% Separate out incorrect should-go ipsis:
%find the cols for no-stim contras that had incorrects (-400 for SC2)
incorrIpsiNoStim = [];
numelColsIpsiNoStim = numel(ipsiNoStimTrajMat(1,:));
incorrIpsiNoStimCol = [];
for i = 1:length(ipsiNoStimTrajMat)                            %for ea row
    for j = 1:numelColsIpsiNoStim                              %for ea col
        if ipsiNoStimTrajMat(i,j) == incorrIpsiT_hold
            incorrIpsiNoStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
incorrIpsiNoStim1 = incorrIpsiNoStimCol > 1;
incorrIpsiNoStim = incorrIpsiNoStimCol(incorrIpsiNoStim1);     %these are the cols of incorr contra trials

%find the cols for stim contras that had incorrects (-400 for SC2)
incorrIpsiStim = [];
numelColsIpsiStim = numel(ipsiStimTrajMat(1,:));
incorrIpsiStimCol = [];
for i = 1:length(ipsiStimTrajMat)                            %for ea row
    for j = 1:numelColsIpsiStim                              %for ea col
        if ipsiStimTrajMat(i,j) == incorrIpsiT_hold
            incorrIpsiStimCol(j) = j;                        %tell me the col where incorr t-hold is reached
        end
    end
end
incorrIpsiStim1 = incorrIpsiStimCol > 1;
incorrIpsiStim = incorrIpsiStimCol(incorrIpsiStim1);         %these are the cols of incorr contra trials

%% Plot the should-move contra trials
figure; hold on;
% plot(contraNoStimTrajMat(:,1:20), timeVec, 'k', 'LineWidth',2);
% lastRow = end-1;

% subplot(2,2,1);
% plot(contraNoStimTrajMat(1:11, 8:10), contraNoStimRelTime(1:11, 8:10), 'bl', 'LineWidth',2);%for blocks
plot(contraNoStimTrajMat(1:27, 16:17), contraNoStimRelTime(1:27, 16:17), 'r', 'LineWidth',2);%for rando

xlabel('position', 'FontSize', 28); ylabel('time (ms) from trial start', 'FontSize', 28);
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:200:1000], 'FontSize', 28);
        ax = gca; 
        ax.FontSize = 28;

hold on;
% set(gca, 'xTickLabel', 'Fontsize', 15, 'yTickLabel', 'FontSize', 15); 
%     mouseInfo = strcat('trials 1:20; should-contra NOstim', ',', mouseName, ',', num2str(numSession));
%     title(mouseInfo);
    
% subplot(2,2,2);
% % plot(contraStimTrajMat(:,[1:8 11:20]), contraStimRelTime(1:478,[1:8 11:20]), 'bl', 'LineWidth',1);
% plot(contraStimTrajMat(:,1:20), contraStimRelTime(1:end-1,1:20), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 1:20;should-contra +stim');
%     title(mouseInfo);

% subplot(2,2,3);
% plot(contraNoStimTrajMat(:,21:40), contraNoStimRelTime(1:end,21:40), 'k', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 21:40;should-contra NOstim');
%     title(mouseInfo);
% subplot (2,2,4)
% plot(contraStimTrajMat(:,21:40), contraStimRelTime(1:end,21:40), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trials 21:40;should-contra +stim');
%     title(mouseInfo);
% 
% %Plot the corrects only (top panels, new fig)
% figure; hold on;
% subplot(2,2,1);
% plot(contraNoStimTrajMat(:,corrContraNoStim), contraNoStimRelTime(1:end,corrContraNoStim), 'k', 'LineWidth',1);
% xlabel('position'); ylabel('time (ms) from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('corrects; should-contra NOstim', ',', mouseName, ',', num2str(numSession));
%     title(mouseInfo);
% subplot(2,2,2);
% plot(contraStimTrajMat(:,corrContraStim), contraStimRelTime(1:end,corrContraStim), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('corrects;should-contra +stim');
%     title(mouseInfo);
% 
% %Plot the incorrects only (bottom panels, new fig)
% subplot(2,2,3);
% plot(contraNoStimTrajMat(:,incorrContraNoStim), contraNoStimRelTime(1:end,incorrContraNoStim), 'k', 'LineWidth',1);
% xlabel('position'); ylabel('time (ms) from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('incorrects; should-contra NOstim', ',', mouseName, ',', num2str(numSession));
%     title(mouseInfo);
% subplot(2,2,4);
% plot(contraStimTrajMat(:,incorrContraStim), contraStimRelTime(1:end,incorrContraStim), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('incorrects;should-contra +stim');
%     title(mouseInfo);
% 

%% Plot the should-move ipsi trials
% figure; hold on;
% subplot(2,2,1); 
plot(ipsiNoStimTrajMat(1:45,16:17), ipsiNoStimRelTime(1:45,16:17), 'r', 'LineWidth',2); %for rando
% plot(ipsiNoStimTrajMat(1:12,48:49), ipsiNoStimRelTime(1:12,48:49), 'bl', 'LineWidth',2); %for blocks
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 1:20;should-ipsi NOstim', ',', mouseName, ',', num2str(numSession));
%     title(mouseInfo);
% subplot(2,2,2); 
% plot(ipsiStimTrajMat(:,1:20), ipsiStimRelTime(1:end-1,1:20), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 1:20;should-ipsi +stim');
%     title(mouseInfo);
% subplot(2,2,3); 
% plot(ipsiNoStimTrajMat(:,21:40), ipsiNoStimRelTime(1:end-1,21:40), 'k', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 21:40;should-ipsi NOstim');
%     title(mouseInfo);
% subplot(2,2,4); 
% plot(ipsiStimTrajMat(:,21:40), ipsiStimRelTime(1:end-1,21:40), 'bl', 'LineWidth',1);
% xlabel('position'); ylabel('ms from trial start');
% set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
%     mouseInfo = strcat('trial 21:40; should-ipsi +stim');
%     title(mouseInfo);
%     
%Plot the corrects only (top panels, new fig)
figure; hold on;
subplot(2,2,1);
plot(ipsiNoStimTrajMat(:,corrIpsiNoStim), ipsiNoStimRelTime(1:end,corrIpsiNoStim), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('corrects; should-ipsi NOstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,2);
plot(ipsiStimTrajMat(:,corrIpsiStim), ipsiStimRelTime(1:end,corrIpsiStim), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('ms from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('corrects;should-ipsi +stim');
    title(mouseInfo);
%Plot the inorrects (bottom panels)
subplot(2,2,3);
plot(ipsiNoStimTrajMat(:,incorrIpsiNoStim), ipsiNoStimRelTime(1:end,incorrIpsiNoStim), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('incorrects; should-ipsi NOstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,4);
plot(ipsiStimTrajMat(:,incorrIpsiStim), ipsiStimRelTime(1:end,incorrIpsiStim), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('ms from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('incorrects;should-ipsi +stim');
    title(mouseInfo);

   
%% Plot all avg traces:
% contra/ipsi all mean trajectories:
contraStimTrajMean = nanmean(contraStimTrajMat(:,1:end),2);
ipsiStimTrajMean = nanmean(ipsiStimTrajMat(:,1:end),2);
contraNoStimTrajMean = nanmean(contraNoStimTrajMat(:,1:end),2);
ipsiNoStimTrajMean = nanmean(ipsiNoStimTrajMat(:,1:end),2);

%is there a significant diff bt mean trajectories?
% [pDist_trajContra,hDist_trajContra,statsDist_contra] = ranksum(contraStimTrajMean(1:20,:), contraNoStimTrajMean(1:20,:)); %these include correct & incorrect trials
% [pDist_trajIpsi,hDist_trajIpsi,statsDist_ipsi] = ranksum(ipsiStimTrajMean(1:20,:), ipsiNoStimTrajMean(1:20,:));           %these include correct & incorrect trials
[pDist_trajContra,hDist_trajContra,statsDist_contra] = ranksum(contraStimTrajMean, contraNoStimTrajMean); %these include correct & incorrect trials
[pDist_trajIpsi,hDist_trajIpsi,statsDist_ipsi] = ranksum(ipsiStimTrajMean, ipsiNoStimTrajMean);           %these include correct & incorrect trials

% time averages:
% contraNoStimRelTimeMean = nanmean(relTimeMat(:,contraNoStim),2);
% contraStimRelTimeMean = nanmean(relTimeMat(:,contraStim),2);
% ipsiNoStimRelTimeMean = nanmean(relTimeMat(:,ipsiNoStim),2);
% ipsiStimRelTimeMean = nanmean(relTimeMat(:,ipsiStim),2);

% Plot the should-move contra trials' average trajectory
figure; hold on;
subplot(2,2,1);
plot(contraNoStimTrajMean, contraNoStimRelTime(1:end,10), 'k', 'LineWidth',1);
% plot(contraNoStimTrajMean, contraNoStimRelTimeMean, 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg trls (all); should-contra NOstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,2);
% plot(contraStimTrajMat(:,[1:8 11:20]), contraStimRelTime(1:478,[1:8 11:20]), 'bl', 'LineWidth',1);
plot(contraStimTrajMean, contraStimRelTime(1:end,9), 'bl', 'LineWidth',1);
% plot(contraStimTrajMean, contraStimRelTimeMean, 'bl', 'LineWidth',1);
xlabel('position'); ylabel('ms from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg trls (all); should-contra +stim');
    title(mouseInfo);

% Plot the should-move ipsi trials' mean
figure; hold on;
subplot(2,2,1); 
plot(ipsiNoStimTrajMean, ipsiNoStimRelTime(1:end,5), 'k', 'LineWidth',1);
xlabel('position'); ylabel('ms from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg trls (all);should-ipsi NOstim', ',', mouseName, ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,2); 
plot(ipsiStimTrajMean, ipsiStimRelTime(1:end,5), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('ms from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg trls (all);should-ipsi +stim');
    title(mouseInfo);

%% Plot the correct v. incorrect means for contra & ipsi:

corrContraStimTrajMean = nanmean(contraStimTrajMat(:,corrContraStim),2); 
incorrContraStimTrajMean = nanmean(contraStimTrajMat(:,incorrContraStim),2);
corrIpsiStimTrajMean = nanmean(ipsiStimTrajMat(:,corrIpsiStim),2);
incorrIpsiStimTrajMean = nanmean(ipsiStimTrajMat(:,incorrIpsiStim),2);

corrContraNoStimTrajMean = nanmean(contraNoStimTrajMat(:,corrContraNoStim),2);
incorrContraNoStimTrajMean = nanmean(contraNoStimTrajMat(:,incorrContraNoStim),2);
corrIpsiNoStimTrajMean = nanmean(ipsiNoStimTrajMat(:,corrIpsiNoStim),2);
incorrIpsiNoStimTrajMean = nanmean(ipsiNoStimTrajMat(:,incorrIpsiNoStim),2);

%Plot the contra trials' correct/incorrect mean trajectories:
figure; hold on;
subplot(2,2,1);
plot(corrContraNoStimTrajMean, contraNoStimRelTime(1:end,corrContraNoStim(1)), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg correct; should-contra NOstim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,2);
plot(corrContraStimTrajMean, contraStimRelTime(1:end,corrContraStim(1)), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg correct; should-contra +stim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,3);
plot(incorrContraNoStimTrajMean, contraNoStimRelTime(1:end,incorrContraNoStim(1)), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg incorr; should-contra NOstim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,4);
plot(incorrContraStimTrajMean, contraStimRelTime(1:end,incorrContraStim(1)), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg incorr; should-contra +stim', ',', num2str(numSession));
    title(mouseInfo);

%Plot the ipsi trials' correct/incorrect mean trajectories:
figure; hold on;
subplot(2,2,1);
plot(corrIpsiNoStimTrajMean, ipsiNoStimRelTime(1:end,corrIpsiNoStim(1)), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg correct; should-Ipsi NOstim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,2);
plot(corrIpsiStimTrajMean, ipsiStimRelTime(1:end,corrIpsiStim(1)), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg correct; should-Ipsi +stim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,3);
plot(incorrIpsiNoStimTrajMean, ipsiNoStimRelTime(1:end,incorrIpsiNoStim(1)), 'k', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg incorr; should-Ipsi NOstim', ',', num2str(numSession));
    title(mouseInfo);
subplot(2,2,4);
plot(incorrIpsiStimTrajMean, ipsiStimRelTime(1:end,incorrIpsiStim(1)), 'bl', 'LineWidth',1);
xlabel('position'); ylabel('time (ms) from trial start');
set(gca, 'Ylim', [0 1000], 'Xlim', [-450 450], 'Ytick', [0:100:1000]);
    mouseInfo = strcat('avg incorr; should-Ipsi +stim', ',', num2str(numSession));
    title(mouseInfo);

