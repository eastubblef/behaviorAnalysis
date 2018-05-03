
%% This script is for generating a heat map (imagesc) of binned rx times' %corrects 

% Calls files from a dir (either for randomized sessions or blocked)
% Bin sizes are currently hard-coded (in the below function that is called)
% 11.2017  This script calls rxSortMat4perform2 or rxSortMat4perform2_perc for forming the matrix

% Updated 12.08.17; Calls these 2 functions:
%   - rxSortMat4perform3.m (updated 12.16.18 to find % priors that were incorr/total incorrect per bin
%   - rxSortOneBack.m for plotting % of trials that were incorrect bc mouse moved in prev rew. dir  
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
  
  [rxSortNHistory] = rxSortOneBack(rxMat);
%   [rxSort] = rxSort4Mat(rxMat);  %call this m file to sort out rx times
%    rxSort = rxSort';              %row 1 = trial num; row 2 = sorted rx times; row 3 = correct L/R v incorrect

   rxSort = rxSortNHistory';      %row 4 = incorrect trials in which priors were rewarded to that dir 
%    rxSort = rxSort';              %row 1 = trial num; row 2 = sorted rx times; row 3 = correct L/R v incorrect
   
%    [Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000] = rxSortMat4perform2(rxSort);
%    [Num100, Num200, Num400, Num600, Num1000, Num1000to3000] = rxSortMat4perform2(rxSort);
%    [Num100, numelCorVec100, Num200, numelCorVec200, Num400, numelCorVec400, Num600, numelCorVec600, Num1000, numelCorVec1000, Num1000to3000, numelCorVec1100] = rxSortMat4perform2(rxSort)
%    [Num100, numelCorVec100, numIncorrPriors100, Num200, numelCorVec200,numIncorrPriors200, Num400, numelCorVec400,numIncorrPriors400, Num600, numelCorVec600,numIncorrPriors600, Num1000, numelCorVec1000,numIncorrPriors1000, Num1000to3000, numelCorVec1100, numIncorrPriors3000] = rxSortMat4perform2(rxSort)
%    [Num100, numelCorVec100, numelIncorVec100, numIncorrPriors100, Num200, numelCorVec200, numelIncorVec200, numIncorrPriors200, Num400, numelCorVec400, numelIncorVec400, numIncorrPriors400,...
%          Num600, numelCorVec600, numelIncorVec600, numIncorrPriors600, Num1000, numelCorVec1000, numelIncorVec1000, numIncorrPriors1000, Num1000to3000, numelCorVec1100, numelIncorVec1100, numIncorrPriors3000] = rxSortMat4perform3(rxSort)        [Num100, numelCorVec100, numelIncorVec100, numIncorrPriors100, percIncorrPriors100ofIncorr100, Num200, numelCorVec200, numelIncorVec200, numIncorrPriors200, percIncorrPriors200ofIncorr200, Num400, numelCorVec400, numelIncorVec400, numIncorrPriors400, percIncorrPriors400ofIncorr400,...

[Num100, numelCorVec100, numelIncorVec100, numIncorrPriors100, percIncorrPriors100ofIncorr100, Num200, numelCorVec200, numelIncorVec200, numIncorrPriors200, percIncorrPriors200ofIncorr200, Num400, numelCorVec400, numelIncorVec400, numIncorrPriors400, percIncorrPriors400ofIncorr400,...
         Num600, numelCorVec600, numelIncorVec600, numIncorrPriors600, percIncorrPriors600ofIncorr600, Num1000, numelCorVec1000, numelIncorVec1000, numIncorrPriors1000, percIncorrPriors1000ofIncorr1000, Num1000to3000, numelCorVec1100, numelIncorVec1100, numIncorrPriors3000, percIncorrPriors3000ofIncorr3000] = rxSortMat4perform3(rxSort)

%   cat = horzcat(Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000);
  cat = horzcat(Num100, Num200, Num400, Num600, Num1000, Num1000to3000);                                                                         %total number of trials per rxn time bin
  catNumTrials = horzcat(numelCorVec100, numelCorVec200, numelCorVec400, numelCorVec600, numelCorVec1000, numelCorVec1100);                      %number of correct trials per rxn time bin
  
  catNumIncorVec = horzcat(numelIncorVec100, numelIncorVec200, numelIncorVec400, numelIncorVec600, numelIncorVec1000, numelIncorVec1100);        %number of incorrect trials per rxn time bin
  catPriors = horzcat(numIncorrPriors100, numIncorrPriors200, numIncorrPriors400, numIncorrPriors600, numIncorrPriors1000, numIncorrPriors3000); %number of trials that were incorrect bc prev. trial was incorrect
  catPriorsPercIncorr = horzcat(percIncorrPriors100ofIncorr100, percIncorrPriors200ofIncorr200, percIncorrPriors400ofIncorr400, percIncorrPriors600ofIncorr600, percIncorrPriors1000ofIncorr1000, percIncorrPriors3000ofIncorr3000);  
 
  rxPerformMat(k,:) = cat;               %out of all trials of these rx times
  rxNumTrialsMat(k,:) = catNumTrials;
  
  rxPriors(k,:) = catPriors;             %out of incorr trials of these rx times
  rxPriorsPercIncorr(k,:) = catPriorsPercIncorr;   %out of incorr trials of these rx times
  rxCatNumIncorVec(k,:) = catNumIncorVec;%number of incorr trials

  clear rxSort;
  clear cat;
  clear catNumTrials;
  
end

%% Josh additionally suggested number of trials & % trials of ea. perform/rxn time bin

% meanRxNumTrialsHist = mean(rxNumTrialsMat, 1);   %mean down ea. col (rx time bin)
numRxBins = 1:length(rxNumTrialsMat(1, :));
sumRxNumTrials = [];
for i = 1:numel(numRxBins)
    sumRxNumTrials(:,i) = sum(rxNumTrialsMat(:,i));
end

%calculate % trials/total for yaxis = %Trials:
sumRxNumTrials_all = sum(sumRxNumTrials);  %total num trials
for l = 1:length(sumRxNumTrials)
    sumRxNumPerc(l) = (sumRxNumTrials(l)/sumRxNumTrials_all) * 100;
end

% To calculate the following, see RxAvg_Performance2.m & rxSortMat4perform2.m - New
%calculate % trials/total per bin in which trial was missed & they went same dir as previously rew. trial
sumRxNumPriors = sum(rxPriors);  
for m = 1:length(sumRxNumTrials)
    sumRxNumPriorsPerc(m) = (sumRxNumPriors(m)/sumRxNumTrials_all) * 100; % percent of total trials mouse moved in prev. rew dir
end

%calculate % trials/incorrects per bin in which trial was missed & they went same dir as previously rew. trial
sumCatNumIncorVec = sum(rxCatNumIncorVec);
sumCatNumIncorVec_all = sum(sumCatNumIncorVec);
for m = 1:length(sumRxNumTrials)
    NumPriorsPercOfIncor(m) = (sumRxNumPriors(m)/sumCatNumIncorVec_all) * 100; % percent of incorrect trials mouse moved in prev. rew dir
end

%% Create bins for rx times and insert % corrects

edges = [100 200 400 600 1000 3000];
figure(1); hold on;
barSize = 250;
% barSize = 150;

for i = 1:length(edges)
    h(i) = bar(edges(i), sumRxNumTrials(i), 'BarWidth', barSize);
%     get(h(i))
%     set(h, 'BarWidth', 250, 'k'); %error
end
set(h, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('# Trials');
set(gca, 'Xlim', [0 3000]) %, 'Ylim', [0 60]);
set(gca, 'Xtick', [0:500:5000]);     
hold off;

barSize = 250;
% barSize = 150;

figure(2); hold on;  %plot y-axis as % Trials
for l = 1:length(edges) 
    i(l) = bar(edges(l), sumRxNumPerc(l), 'BarWidth', barSize);
%     set(h(l), 'BarWidth', barSize);
end
set(i, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('% Trials');
set(gca, 'Xlim', [0 3000], 'Ylim', [0 60]);
set(gca, 'Xtick', [0:500:5000]);     barSize = 250;

hold off;

%for priors:

figure(3); hold on;  %plot y-axis as % Trials
% for l = 1:length(edges) 
%     h(l) = bar(edges(l), sumRxNumPriorsPerc(l));
%     set(h(l), 'BarWidth', barSize);
% end
for l = 1:length(edges) 
    j(l) = bar(edges(l), NumPriorsPercOfIncor(l), 'BarWidth', barSize);
    set(h(l), 'BarWidth', barSize);
end
set(j, 'FaceColor', 'k');
xlabel('Rxn times (ms)'); ylabel('% Incorrect trials');
set(gca, 'Xlim', [0 3000]), 
set(gca, 'Xtick', [0:500:5000]);    

hold off;

figure(4); hold on;
h_image = imagesc(rxPerformMat);
colorbar;
ylabel('Randomized sessions'); xlabel('Rxn times(edges)'); 

%% For additional averaging down the columns of reaction times:

%% For additional averaging down the columns of reaction times:

rxPerformMatAvg = nanmean(rxPerformMat, 1);
rxPerformMatAvg = round(rxPerformMatAvg);
figure; hold on;
h_imageAvg_blks = imagesc(rxPerformMatAvg(1:6)); %blocks % corrects
%From Josh:
[cmap] = TNC_CreateRBColormap(1024,'gp');
colormap(cmap);
c_blcks = colorbar;
cmin = 0; cmax = 100; caxis([cmin cmax]);

figure; hold on;
h_hardCode = imagesc([64 64 87 95 89 86]); %randomized % corrects

[cmap] = TNC_CreateRBColormap(1024,'gp');
colormap(cmap);
c_rand = colorbar;
cmin = 0; cmax = 100; caxis([cmin cmax]);
% % To update the colormap being applies to a figure call:
% colormap(cmap);
% hh = imagesc(j,'CDataMapping','scaled')

% Break down ea. performance column color to apply to the histograms:
% h_imageAvg100 = imagesc(rxPerformMatAvg(1));
% set(h_image.YData)
% h_imageAvg100_200 = imagesc(rxPerformMatAvg(2));
% h_imageAvg200_400 = imagesc(rxPerformMatAvg(3));
% h_imageAvg400_600 = imagesc(rxPerformMatAvg(4));
% h_imageAvg600_1000 = imagesc(rxPerformMatAvg(5));
% h_imageAvg1000_3000 = imagesc(rxPerformMatAvg(6));








