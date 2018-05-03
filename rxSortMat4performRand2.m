% function  [Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000] = rxSortMat4perform2(rxSort)
function   [Num100, Num200, Num400, Num600, Num1000, Num1000to3000] = rxSortMat4perform2Rand(rxSort)


%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% edges = [1:5:3000];                  %gives 600 bins for Vgatfive 12.15.16; 23 vals in first bin: numel(find(rxSort(:,2) < 5));
% edges = [5:10:3000];                 %gives 300 bins for Vgatfive 12.15.16; 23 vals
edges = [0:100:3000];                  %gives 27 bins for Vgatfive 12.15.16
% numBins = numel(edges);
mids = conv2(edges, [1 1], 'valid')/2; %gives 26 bins
numBins = numel(edges);

%Get the spread of incorrects/corrects within bins by indexing back into rxsort 
mu = rxSort(3,:);
rxSort1 = rxSort(2,:);

binwidth = edges(2) - edges(1);
binwidth100 = binwidth;
hundredInds = find(rxSort1 <= binwidth100);  
rxSort100 = rxSort1(hundredInds);
corVec100 = mu(hundredInds);
yes100 = find(corVec100 ~= 0); 
numelCorVec100 = numel(corVec100);
perc100Cor = numel(yes100)/numelCorVec100;   %percent corrects for fastest bin
Num100 = perc100Cor * 100;
% rxSortFastest = rxSort1(1:binwidthFastest);              %rx times for fastest bin

% binwidthFast = mids(4) - mids(1);
binwidth200 = edges(3) 
twoHunInds = find(rxSort1 > edges(2) & rxSort1 <= edges(3));  
rxSort200 = rxSort1(twoHunInds);                          %rx times for fast bin (200-600 ms)
corVec200 = mu(rxSort1 > edges(2) & rxSort1 <= edges(3));
yes200 = find(corVec200 ~=0);
numelCorVec200 = numel(corVec200);
perc200Cor = numel(yes200)/numelCorVec200;
Num200 = perc200Cor * 100;                             %percent corrects for fast bin

% binwidth300 = edges(4);
% threeHunInds = find(rxSort1 > edges(3) & rxSort1 <= edges(4));  
% rxSort300 = rxSort1(threeHunInds);                           %rx times for mid bin (600-1000 ms)
% corVec300 = mu(rxSort1 > edges(3) & rxSort1 <= edges(4));  
% yes300 = find(corVec300 ~= 0);
% numelCorVec300 = numel(corVec300);
% perc300Cor = numel(yes300)/numelCorVec300;
% Num300 = perc300Cor * 100;                              %percent corrects for mid bin

binwidth400 = edges(5) %taylor to 200-400ms
fourHunInds = find(rxSort1 > edges(3) & rxSort1 <= edges(5)); 
rxSort400 = rxSort1(fourHunInds);                         %rx times for the mid-range bin (600-1000 ms)
corVec400 = mu(rxSort1 > edges(3) & rxSort1 <= edges(5));
yes400 = find(corVec400 ~= 0);
numelCorVec400 = numel(corVec400);
perc400Cor = numel(yes400)/numelCorVec400;
Num400 = perc400Cor * 100;                            %percent corrects for slow bin

% binwidth500 = edges(6)
% fiveHunInds = find(rxSort1 > edges(5) & rxSort1 <= edges(6)); 
% rxSort500 = rxSort1(fiveHunInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec500 = mu(rxSort1 > edges(5) & rxSort1 <= edges(6)); ;
% yes500 = find(corVec500 ~= 0);
% numelCorVec500 = numel(corVec500);
% perc500Cor = numel(yes500)/numelCorVec500;
% Num500 = perc500Cor * 100;                            %percent corrects for slow bin

binwidth600 = edges(7)  %tailor to 400-600
sixHunInds = find(rxSort1 > edges(5) & rxSort1 <= edges(7)); 
rxSort600 = rxSort1(sixHunInds);                         %rx times for the mid-range bin (600-1000 ms)
corVec600 = mu(rxSort1 > edges(5) & rxSort1 <= edges(7)); 
yes600 = find(corVec600 ~= 0);
numelCorVec600 = numel(corVec600);
perc600Cor = numel(yes600)/numelCorVec600;
Num600 = perc600Cor * 100;                            %percent corrects for slow bin

% binwidth700 = edges(8)
% sevHunInds = find(rxSort1 > edges(7) & rxSort1 <= edges(8)); 
% rxSort700 = rxSort1(sevHunInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec700 = mu(rxSort1 > edges(7) & rxSort1 <= edges(8)); 
% yes700 = find(corVec700 ~= 0);
% numelCorVec700 = numel(corVec700);
% perc700Cor = numel(yes700)/numelCorVec700;
% Num700 = perc700Cor * 100;                            %percent corrects for slow bin
% 
% binwidth800 = edges(9)
% eightHunInds = find(rxSort1 > edges(8) & rxSort1 <= edges(9)); 
% rxSort800 = rxSort1(eightHunInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec800 = mu(rxSort1 > edges(8) & rxSort1 <= edges(9)); 
% yes800 = find(corVec800 ~= 0);
% numelCorVec800 = numel(corVec800);
% perc800Cor = numel(yes800)/numelCorVec800;
% Num800 = perc800Cor * 100;                            %percent corrects for slow bin
% 
% binwidth900 = edges(10)
% nineHunInds = find(rxSort1 > edges(9) & rxSort1 <= edges(10)); 
% rxSort900 = rxSort1(nineHunInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec900 = mu(rxSort1 > edges(9) & rxSort1 <= edges(10)); 
% yes900 = find(corVec900 ~= 0);
% numelCorVec900 = numel(corVec900);
% perc900Cor = numel(yes900)/numelCorVec900;
% Num900 = perc900Cor * 100;                            %percent corrects for slow bin

binwidth1000 = edges(11);  %tailor to 600-1000ms
thouInds = find(rxSort1 > edges(7) & rxSort1 <= edges(11)); 
rxSort1000 = rxSort1(thouInds);                         %rx times for the mid-range bin (600-1000 ms)
corVec1000 = mu(rxSort1 > edges(7) & rxSort1 <= edges(11)); 
yes1000 = find(corVec1000 ~= 0);
numelCorVec1000 = numel(corVec1000);
perc1000Cor = numel(yes1000)/numelCorVec1000;
Num1000 = perc1000Cor * 100;                            %percent corrects for slow bin

binwidthLargerThan1000 = edges(11) %tailor to 1000-3000
thouoneInds = find(rxSort1 > edges(11) & rxSort1 <= edges(31));  
rxSort1100 = rxSort1(thouoneInds);                         %rx times for the mid-range bin (600-1000 ms)
corVec1100 = mu(rxSort1 > edges(11) & rxSort1 <= edges(31)); 
yes1100 = find(corVec1100 ~= 0);
numelCorVec1100 = numel(corVec1100);
perc1100Cor = numel(yes1100)/numelCorVec1100;
Num1000to3000 = perc1100Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1200 = edges(13)
% thoutwoInds = find(rxSort1 > edges(12) & rxSort1 <= edges(13)); 
% rxSort1200 = rxSort1(thoutwoInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1200 = mu(rxSort1 > edges(12) & rxSort1 <= edges(13)); 
% yes1200 = find(corVec1200 ~= 0);
% numelCorVec1200 = numel(corVec1200);
% perc1200Cor = numel(yes1200)/numelCorVec1200;
% Num1200 = perc1200Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1300 = edges(14)
% thouthreeInds = find(rxSort1 > edges(13) & rxSort1 <= edges(14)); 
% rxSort1300 = rxSort1(thouthreeInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1300 = mu(rxSort1 > edges(13) & rxSort1 <= edges(14)); 
% yes1300 = find(corVec1300 ~= 0);
% numelCorVec1300 = numel(corVec1300);
% perc1300Cor = numel(yes1300)/numelCorVec1300;
% Num1300 = perc1300Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1400 = edges(15)
% thoufourInds = find(rxSort1 > edges(14) & rxSort1 <= edges(15)); 
% rxSort1400 = rxSort1(thoufourInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1400 = mu(rxSort1 > edges(14) & rxSort1 <= edges(15)); 
% yes1400 = find(corVec1400 ~= 0);
% numelCorVec1400 = numel(corVec1400);
% perc1400Cor = numel(yes1400)/numelCorVec1400;
% Num1400 = perc1400Cor * 100;                            %percent corrects for slow bin

% binwidth1500 = edges(16)  %tailor to 1500-2000
% thoufiveInds = find(rxSort1 > edges(16) & rxSort1 <= edges(21)); 
% rxSort1500 = rxSort1(thoufiveInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1500 = mu(rxSort1 > edges(16) & rxSort1 <= edges(21)); 
% yes1500 = find(corVec1500 ~= 0);
% numelCorVec1500 = numel(corVec1500);
% perc1500Cor = numel(yes1500)/numelCorVec1500;
% Num1500 = perc1500Cor * 100;                            %percent corrects for slow bin

% binwidth1600 = edges(17)
% thousixInds = find(rxSort1 > edges(16) & rxSort1 <= edges(17)); 
% rxSort1600 = rxSort1(thousixInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1600 = mu(rxSort1 > edges(16) & rxSort1 <= edges(17)); 
% yes1600 = find(corVec1600 ~= 0);
% numelCorVec1600 = numel(corVec1600);
% perc1600Cor = numel(yes1600)/numelCorVec1600;
% Num1600 = perc1600Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1700 = edges(18)
% thousevInds = find(rxSort1 > edges(17) & rxSort1 <= edges(18)); 
% rxSort1700 = rxSort1(thousevInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1700 = mu(rxSort1 > edges(17) & rxSort1 <= edges(18)); 
% yes1700 = find(corVec1700 ~= 0);
% numelCorVec1700 = numel(corVec1700);
% perc1700Cor = numel(yes1700)/numelCorVec1700;
% Num1700 = perc1700Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1800 = edges(19)
% thoueightInds = find(rxSort1 > edges(18) & rxSort1 <= edges(19)); 
% rxSort1800 = rxSort1(thoueightInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1800 = mu(rxSort1 > edges(18) & rxSort1 <= edges(19)); 
% yes1800 = find(corVec1800 ~= 0);
% numelCorVec1800 = numel(corVec1800);
% perc1800Cor = numel(yes1800)/numelCorVec1800;
% Num1800 = perc1800Cor * 100;                            %percent corrects for slow bin
% 
% binwidth1900 = edges(20)
% thounineInds = find(rxSort1 > edges(19) & rxSort1 <= edges(20)); 
% rxSort1900 = rxSort1(thounineInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec1900 = mu(rxSort1 > edges(19) & rxSort1 <= edges(20)); 
% yes1900 = find(corVec1900 ~= 0);
% numelCorVec1900 = numel(corVec1900);
% perc1900Cor = numel(yes1900)/numelCorVec1900;
% Num1900 = perc1900Cor * 100;                            %percent corrects for slow bin

% binwidth2000 = edges(16:21)
% twothouInds = find(rxSort1 > edges(16) & rxSort1 <= edges(21)); 
% rxSort2000 = rxSort1(twothouInds);                         %rx times for the mid-range bin (600-1000 ms)
% corVec2000 = mu(rxSort1 > edges(16) & rxSort1 <= edges(21)); 
% yes2000 = find(corVec2000 ~= 0);
% numelCorVec2000 = numel(corVec2000);
% perc2000Cor = numel(yes2000)/numelCorVec2000;
% Num2000 = perc2000Cor * 100;                            %percent corrects for slow bin



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

