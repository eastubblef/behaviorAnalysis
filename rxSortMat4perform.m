function [NumFastest, NumFast, NumMid, NumSlow] = rxSortMat4perform(rxSort)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% edges = [1:5:3000];                  %gives 600 bins for Vgatfive 12.15.16; 23 vals in first bin: numel(find(rxSort(:,2) < 5));
% edges = [5:10:3000];                 %gives 300 bins for Vgatfive 12.15.16; 23 vals
edges = [0:200:5200];                  %gives 27 bins for Vgatfive 12.15.16
% numBins = numel(edges);
mids = conv2(edges, [1 1], 'valid')/2; %gives 26 bins
numBins = numel(edges);

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
mu = rxSort(3,:);
rxSort1 = rxSort(2,:);

binwidth = edges(2) - edges(1);
binwidthFastest = edges(2) - edges(1);

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

%% create a matrix of rx times in 100 ms bins (eventually) for ea. session and put % corrects in each 


%% Do the plotting
% % histRx = histc(rxSort(:,2), edges);                    % all rx times
% histRx1 = histc(rxSortFastest, edges);                   % 0-200 ms bin
% histRx2 = histc(rxSortFast, edges);                      % 200-600 ms bin
% histRx3 = histc(rxSortMid, edges);                       % 600-1000 ms bin
% histRx4 = histc(rxSortSlow, edges);                      % >1000 ms bin
% 
% figure; hold on;
% h4 = bar(edges, histRx4, 'k');
% h3 = bar(edges, histRx3, 'g'); 
% h2 = bar(edges, histRx2, 'b');
% h1 = bar(edges, histRx1, 'r'); set(h1, 'BarWidth', 1);  % plot last for visual purposes

%     xlabel('Rx times (ms), 200ms bins'); ylabel('# Trials');
% %   set(gca, 'Xlim', [min(edges) max(edges)]);
%     set(gca, 'Xlim', [0 5200]);
% %   set(gca, 'Xtick', [min(edges):600 max(edges)]), %'Xticklabel', [min(edges):max(edges)]);
%     set(gca, 'Xtick', [0:500:5200]), 
%     num200 = num2str(round(NumFastest)); num200_600 = num2str(round(NumFast)); num600_1000 = num2str(round(NumMid)); num1000 = num2str(round(NumSlow));
%     legend(num1000,num600_1000,num200_600,num200,'Location','northeast');
%     mouseInfo = strcat('% corrects ', ',', mouseName, ',', num2str(numSessions), ' sessions');
%     title(mouseInfo);
   
end

