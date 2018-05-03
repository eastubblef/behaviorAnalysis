% This script is for extracting the average reaction time per session and
% plotting it with respect to training day.

%% Inputs:

fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbBlocks';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
cd(fpath);
% folder = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
% fullfile = strcat(folder,'/',file);

% numfiles = 70;
rxPerformMat = []
dinfo = dir('*.mat');

for k = 1 : length(dinfo)         %use this for all 70ish sessions
  thismat = dinfo(k).name;
  st = load(thismat);
  rxMat = st.taskbase.rxMat;
  
  sumRx1 = vertcat(rxMat(:,2), rxMat(:,3), rxMat(:,4), rxMat(:,5));
 
  meanRx = nanmean(sumRx1);
end
