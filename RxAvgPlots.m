
%% This script plots the rx times from multiple sessions into one histogram;
% Bin sizes are set for 200 ms = fastest bin; 200-600 ms = fast bin; 600-1000 = mid bin; >1000 = slow bin

%% Inputs:

taskbaseFile1 = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2/mSC2_2017_03_02_14_tb.mat';
taskbaseFile2 = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2/mSC2_2017_03_03_16_tb.mat';
taskbaseFile3 = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2/mSC2_2017_03_06_14_tb.mat';
taskbaseFile4 = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2/mSC2_2017_02_28_15_tb.mat';

mouseName1 = 'SC'; mouseName2 =  num2str(2);  %mouse SC2
% mouseName1 = 'SC'; mouseName2 =  num2str(1);  %mouse SC1
numSessions = 4;

mouseName = strcat(mouseName1, mouseName2);
load(taskbaseFile1);
rxMat1 = taskbase.rxMat;
clear taskbaseFile1;

load(taskbaseFile2);
rxMat2 = taskbase.rxMat;
clear taskbaseFile2;

load(taskbaseFile3);
rxMat3 = taskbase.rxMat;
clear taskbaseFile3;

load(taskbaseFile4);
rxMat4 = taskbase.rxMat;
clear taskbaseFile4;

rxMat = vertcat(rxMat1, rxMat2, rxMat3, rxMat4);

%% Use the matrix of rx times:     col 6: CL = -1; CR = 1; IL = 0; IR = 0;
%col1 = trial number
%col2 = rx times for CL trials
%col3 = rx times for CR trials
%col4 = rx times for IL trials
%col5 = rx times for IR trials

%bin the rx times:
% hist(x,xbins)
rxCondensed = [];
for m = 1:length(rxMat(:,2))
    if isfinite(rxMat(m,2))
        rxCondensed(m,1) = rxMat(m,1);
        rxCondensed(m,2) = rxMat(m,2);
        rxCondensed(m,3) = rxMat(m,6);
    else if isfinite(rxMat(m,3))
            rxCondensed(m,1) = rxMat(m,1);
            rxCondensed(m,2) = rxMat(m,3);
            rxCondensed(m,3) = rxMat(m,6);
        else if isfinite(rxMat(m,4))
                rxCondensed(m,1) = rxMat(m,1);
                rxCondensed(m,2) = rxMat(m,4);
                rxCondensed(m,3) = rxMat(m,6);
            else if isfinite(rxMat(m,5))
                    rxCondensed(m,1) = rxMat(m,1);
                    rxCondensed(m,2) = rxMat(m,5);
                    rxCondensed(m,3) = rxMat(m,6);
                else rxCondensed(m,1) = rxMat(m,1);
                     rxCondensed(m,2) = NaN;
                     rxCondensed(m,3) = rxMat(m,6);
                end
            end
        end
    end
end

% Sorted Rx times & drag trial info with ea:
rx = rxCondensed(:,2);
corrVals = rxCondensed(:,3);
[rxSort1, index] = sort(rx, 1, 'ascend');
rxSort = [];
rxSort(:,1) = index;
rxSort(:,2) = rxSort1;
rxSort(:,3) = corrVals(index);

% edges = [1:5:3000];                  %gives 600 bins for Vgatfive 12.15.16; 23 vals in first bin: numel(find(rxSort(:,2) < 5));
% edges = [5:10:3000];                 %gives 300 bins for Vgatfive 12.15.16; 23 vals
edges = [0:200:5200];                  %gives 27 bins for Vgatfive 12.15.16
% numBins = numel(edges);
mids = conv2(edges, [1 1], 'valid')/2; %gives 26 bins
numBins = numel(edges);

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
mu = rxSort(:,3);
binwidth = edges(2) - edges(1);
binwidthFastest = edges(2) - edges(1);
% corVecFastest = [];
% for i = 1:length(rxSort(:,1))
%     if rxSort1(i) < binwidthFastest 
%         corVecFastest(i) = mu(i);
%     end
% end
% corVecFastest = corVecFastest';

fastestInds = find(rxSort1 < binwidthFastest);  
rxSortFastest = rxSort1(fastestInds);
% corVecFastest = mu(rxSort1 > edges(1) & rxSort1< binwidthFastest);
corVecFastest = mu(fastestInds);

yesFastest = find(corVecFastest ~= 0); 
numelCorVecFastest = numel(corVecFastest);
percFastestCor = numel(yesFastest)/numelCorVecFastest;   %percent corrects for fastest bin
NumFastest = percFastestCor * 100;
% rxSortFastest = rxSort1(1:binwidthFastest);              %rx times for fastest bin

% binwidthFast = mids(4) - mids(1);
binwidthFast = edges(4) - edges(1);
fastInds = find(rxSort1 > edges(2) & rxSort1< binwidthFast);  
rxSortFast = rxSort1(fastInds);                          %rx times for fast bin (200-600 ms)
corVecFast = mu(rxSort1 > edges(2) & rxSort1< binwidthFast);
yesFast = find(corVecFast ~=0);
numelCorVecFast = numel(corVecFast);
percFastCor = numel(yesFast)/numelCorVecFast;
NumFast = percFastCor * 100;                             %percent corrects for fast bin

binwidthMid = edges(6);
midInds = find(rxSort1 > edges(4) & rxSort1 < binwidthMid);  
rxSortMid = rxSort1(midInds);                           %rx times for mid bin (600-1000 ms)
corVecMid = mu(rxSort1 > edges(4) & rxSort1 < binwidthMid);  
yesMid = find(corVecMid ~= 0);
numelCorVecMid = numel(corVecMid);
percMidCor = numel(yesMid)/numelCorVecMid;
NumMid = percMidCor * 100;                              %percent corrects for mid bin

slowInds = find(rxSort1 > binwidthMid);
rxSortSlow = rxSort1(slowInds);                         %rx times for the mid-range bin (600-1000 ms)
corVecSlow = mu(rxSort1 > binwidthMid);
yesSlow = find(corVecSlow ~= 0);
numelCorVecSlow = numel(corVecSlow);
percSlowCor = numel(yesSlow)/numelCorVecSlow;
NumSlow = percSlowCor * 100;                            %percent corrects for slow bin

%% Do the plotting
% histRx = histc(rxSort(:,2), edges);                    % all rx times
histRx1 = histc(rxSortFastest, edges);                   % 0-200 ms bin
histRx2 = histc(rxSortFast, edges);                      % 200-600 ms bin
histRx3 = histc(rxSortMid, edges);                       % 600-1000 ms bin
histRx4 = histc(rxSortSlow, edges);                      % >1000 ms bin

figure; hold on;
h4 = bar(edges, histRx4, 'k');
h3 = bar(edges, histRx3, 'g'); 
h2 = bar(edges, histRx2, 'b');
h1 = bar(edges, histRx1, 'r'); set(h1, 'BarWidth', 1);  % plot last for visual purposes

    xlabel('Rx times (ms), 200ms bins'); ylabel('# Trials');
%   set(gca, 'Xlim', [min(edges) max(edges)]);
    set(gca, 'Xlim', [0 5200]);
%   set(gca, 'Xtick', [min(edges):600 max(edges)]), %'Xticklabel', [min(edges):max(edges)]);
    set(gca, 'Xtick', [0:500:5200]), 
    num200 = num2str(round(NumFastest)); num200_600 = num2str(round(NumFast)); num600_1000 = num2str(round(NumMid)); num1000 = num2str(round(NumSlow));
    legend(num1000,num600_1000,num200_600,num200,'Location','northeast');
    mouseInfo = strcat('% corrects ', ',', mouseName, ',', num2str(numSessions), ' sessions');
    title(mouseInfo);
   