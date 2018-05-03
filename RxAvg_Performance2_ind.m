
%% This script is for generating a heat map (imagesc) of binned rx times' %corrects 

% Calls files from a dir (either for randomized sessions or blocked)
% Bin sizes are currently hard-coded (in the below function that is called)
% 11.2017  This script calls rxSortMat4perform2 or rxSortMat4perform2_perc for forming the matrix
%% Inputs:

fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbBlocks';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbBlocks/SC1';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbBlocks/Vgatfour';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbBlocks/Vgatseven';

% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand/Vgatfive';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand/Vgateight';
% fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand/SC2';

cd(fpath);
% folder = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
% fullfile = strcat(folder,'/',file);

% numfiles = 70;
rxPerformMat = [];

dinfo = dir('*.mat');
% for k = 1 : length(dinfo(1:20))          %use this per mouse (20ish sessions) - Vgatfive
%  for k = 1 : length(dinfo(26:39))        %use this per mouse (20ish sessions) - SC2
%  for k = 1 : length(dinfo(40:70))        %use this per mouse (20ish sessions) - Vgateight

%  for k = 1 : length(dinfo(1:54))         %use this per mouse (20ish sessions) - Vgatseven
%  for k = 1 : length(dinfo(55:64))        %use this per mouse (20ish sessions) - Vgatfour
%  for k = 1 : length(dinfo(65:75))        %use this per mouse (20ish sessions) - SC6

%  for k = 1 : length(dinfo(65:75))        %use this per mouse (20ish sessions) - SC1
%  for k = 1 : length(dinfo(76:108))       %use this per mouse (20ish sessions) - SC1


for k = 1 : length(dinfo)                %use this for all 70ish sessions
  thismat = dinfo(k).name;
  st = load(thismat);
  rxMat = st.taskbase.rxMat;
  
  [rxSortNHistory] = rxSortOneBack(rxMat);
%   [rxSort] = rxSort4Mat(rxMat);  %call this m file to sort out rx times
%    rxSort = rxSort';             %row 1 = trial num; row 2 = sorted rx times; row 3 = correct L/R v incorrect

   rxSort = rxSortNHistory';       %row 4 = incorrect trials in which priors were rewarded to that dir 
%    rxSort = rxSort';             %row 1 = trial num; row 2 = sorted rx times; row 3 = correct L/R v incorrect
   
%    [Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000] = rxSortMat4perform2(rxSort);
%    [Num100, Num200, Num400, Num600, Num1000, Num1000to3000] = rxSortMat4perform2(rxSort);

%    [Num100, numelCorVec100, Num200, numelCorVec200, Num400, numelCorVec400, Num600, numelCorVec600, Num1000, numelCorVec1000, Num1000to3000, numelCorVec1100] = rxSortMat4perform2(rxSort)
     [Num100, numelCorVec100, numIncorrPriors100, Num200, numelCorVec200,numIncorrPriors200, Num400, numelCorVec400,numIncorrPriors400, Num600, numelCorVec600,numIncorrPriors600, Num1000, numelCorVec1000,numIncorrPriors1000, Num1000to3000, numelCorVec1100, numIncorrPriors3000] = rxSortMat4perform2(rxSort)

 %   cat = horzcat(Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000);
  cat = horzcat(Num100, Num200, Num400, Num600, Num1000, Num1000to3000);
  catNumTrials = horzcat(numelCorVec100, numelCorVec200, numelCorVec400, numelCorVec600, numelCorVec1000, numelCorVec1100);
  catPriors = horzcat(numIncorrPriors100, numIncorrPriors200, numIncorrPriors400, numIncorrPriors600, numIncorrPriors1000,numIncorrPriors3000);  
  
  rxPerformMat(k,:) = cat;               %out of all trials of these rx times
  rxNumTrialsMat(k,:) = catNumTrials;
  rxPriors(k,:) = catPriors;         %out of all trials of these rx times
  
  clear rxSort;
  clear cat;
  clear catNumTrials;
  
end

%% Josh additionally wants number of trials of ea.perform/rxn time bin

% meanRxNumTrialsHist = mean(rxNumTrialsMat, 1);   %mean down ea. col (rx time bin)
numRxBins = 1:length(rxNumTrialsMat(1, :));
sumRxNumTrials = [];
for i = 1:numel(numRxBins)
    sumRxNumTrials(:,i) = sum(rxNumTrialsMat(:,i));
end

%calculate % trials/total per bin for yaxis = %Trials:
sumRxNumTrials_all = sum(sumRxNumTrials);  %total num trials
for l = 1:length(sumRxNumTrials)
    sumRxNumPerc(l) = (sumRxNumTrials(l)/sumRxNumTrials_all) * 100;
end

%calculate % trials/total per bin in which trial was missed & they went same dir as previously rew. trial
sumRxNumPriors = sum(rxPriors);  %total num trials
for m = 1:length(sumRxNumTrials)
    sumRxNumPriorsPerc(m) = (sumRxNumPriors(m)/sumRxNumTrials_all) * 100;
end


%% Create bins for rx times and insert % corrects

edges = [100 200 400 600 1000 3000];
figure; hold on;
barSize = 250;
for i = 1:length(edges)
    h(i) = bar(edges(i), sumRxNumTrials(i), 'BarWidth', barSize);
%     get(h(i))
%     set(h, 'BarWidth', 250, 'k'); %error
end
set(h, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('# Trials');
set(gca, 'Xlim', [0 3000]), %'Ylim', [0 50]);
set(gca, 'Xtick', [0:500:5000]);     
hold off;

barSize = 250;
figure; hold on;  %plot y-axis as % Trials
for l = 1:length(edges) 
    h(l) = bar(edges(l), sumRxNumPerc(l));
    set(h(l), 'BarWidth', barSize);
end
set(h, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('% Trials');
set(gca, 'Xlim', [0 3000]), 
set(gca, 'Xtick', [0:500:5000]);     
hold off;

%for priors:

figure; hold on;  %plot y-axis as % Trials
for l = 1:length(edges) 
    h(l) = bar(edges(l), sumRxNumPriorsPerc(l));
    set(h(l), 'BarWidth', barSize);
end
set(h, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('% Trials');
set(gca, 'Xlim', [0 3000]), 
set(gca, 'Xtick', [0:500:5000]);     
hold off;

figure; hold on;
h_image = imagesc(rxPerformMat);
colorbar;
ylabel('Ind. blocked sessions'); xlabel('Rxn times(edges)'); 

