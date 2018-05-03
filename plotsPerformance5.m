%% This script plots total percent performance for a behavioral session.
%Just save your excel file as excel 5.0/95 then it will work - Mathworks
%PLOTS ALL MICE TOGETHER for new coupling only!

%Updated 1.17.18 for better mean overlay 

fileDirectory = '/Users/stubblefielde/Desktop/mfiles/behavior/';
% file = 'mouseBehaviorWorks.xls';  %worked 2.17
file = 'mouseBehavior2.xls';        %updated 9.15.17
filename = strcat(fileDirectory, file);

mousename = 'Vgatfour';             %input the relevant rows from xls spreadsheet
% timespan = 136:218;                 
timespan = 256:329;                 %updated 9.15.17
blocks = 1;
fullbar = 1;
newCoupling = 0;                    %drk blue - old coupling

%% Begin new coupling:
mousename1 = 'Vgatfive';
% timespan1 = 335:413;              %full bar
timespan1 = 338:528;                %full bar updated 9.15.17
blocks1 = 0;
fullbar1 = 1;
newCoupling1 = 0;    %drk red

timespan1two = 529:591;             %half bar up to recording day2
blocks1two = 0;
fullbar1two = 0
newCoupling1two = 0; %magenta

mousename2 = 'Vgatseven'; 
% timespan2 = 688:731                 %full bar, old coupling
timespan2 = 800:856                 %full bar, old coupling
blocks2 = 1;
fullbar2 = 1,                       %solid

% timespan2two = 742:824;             %half bar, new coupling
timespan2two = 857:965;               %updated 9.15.17
blocks2two = 1;
fullbar2two = 0;
newCoupling2two = 1; %cyan

mousename3 = 'Vgateight';
% timespan3 = 838:911                 %includes blocks, starting out, then conversion to randomized
timespan3 = 981:1103;                  %update; includes blocks, starting out, then conversion to randomized
blocks3 = 0;
fullbar3 = 0;
newCoupling3 = 1;   %magenta

mousename4 = 'SC1';
% timespan4 = 926:940;
timespan4 = 1116:1222;                %update
blocks4 = 1;
fullbar4 = 0;
newCoupling4 = 1;   %cyan

mousename5 = 'SC2';
% timespan5 = 953:963;                   %953:967 for full blocks timespan
timespan5 = 1278:1425;                   %update
% blocks5 = 1;
blocks5 = 0;
fullbar5 = 0;
newCoupling5 = 1;   %cyan

timespan5two = 1298:1400; %switching to randomized
blocks5two = 0;
fullbar5two = 0;
newCoupling5two = 1;  %magenta

timespan5two = 1298:1400; %switching to randomized
blocks5two = 0;
fullbar5two = 0;
newCoupling5two = 1;  %magenta


mousename6 = 'SC3';
timespan6 = 1457:1503;                   %update
blocks6 = 1;
fullbar6 = 0;
newCoupling6 = 1;   %cyan


%% Make it an xlsRead function: 

[num, txt] = xlsread(filename); 

%% Plot like Fig 1 of Burgess et al., 2016 
% mouse Vgatfive randomized

trainingDay1 = num(timespan1, 5);
percCorrectsTot1 = num(timespan1, 22);

trainingDay1two = num(timespan1two, 5);
percCorrectsTot1two = num(timespan1two, 22);

percCorrectsMat1(:,1) = trainingDay1;
percCorrectsMat1(:,2) = percCorrectsTot1;

percCorrectsMat1two(:,1) = trainingDay1two;
percCorrectsMat1two(:,2) = percCorrectsTot1two;

logMat1 = isfinite(percCorrectsMat1);
for i = 1:length(percCorrectsMat1)
    if logMat1(i,1) == 1 && logMat1(i,2) == 1         %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat1(i,:) = percCorrectsMat1(i,:);
    end
end
logMat1two = isfinite(percCorrectsMat1two);
for i2 = 1:length(percCorrectsMat1two)
    if logMat1two(i2,1) == 1 && logMat1two(i2,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat1two(i2,:) = percCorrectsMat1two(i2,:);
    end
end

vecTrain1 = find(percCorrectsTotMat1(:,1) > 0);
vecPerform1 = find(percCorrectsTotMat1(:,2) > 0);
performMat1(:,1) = percCorrectsTotMat1(vecTrain1,1);
performMat1(:,2) = percCorrectsTotMat1(vecPerform1,2);
performMat1(:,3) = 1:length(performMat1(:,2));           %normalize to the new task's number of training days

vecTrain1two = find(percCorrectsTotMat1two(:,1) > 0);
vecPerform1two = find(percCorrectsTotMat1two(:,2) > 0);
performMat1two(:,1) = percCorrectsTotMat1two(vecTrain1two,1);
performMat1two(:,2) = percCorrectsTotMat1two(vecPerform1two,2);
performMat1two(:,3) = 1:length(performMat1two(:,2));

vec1one = performMat1(:,1); vec1two = performMat1(:,2);
vec1three = performMat1two(:,1); vec1four = performMat1two(:,2);
time2perform1one1 = find(vec1two == 75);
time2perform1two1 = find(vec1four == 75);

time2perform1one = vec1one(time2perform1one1);
time2perform1two = vec1three(time2perform1two1);

% Plot 
figure; hold on;
%     ax = gca;
%     mn = min(performMat1(:,1)); mx = max(performMat1two(:,1));
%     y1 = 75; y2 = 75;
%     l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
%     set(l,'LineStyle','--')

% plot(performMat1(:,1), performMat1(:,2), 'm', 'LineWidth', 2);
p = plot(performMat1two(:,3), performMat1two(:,2), 'y', 'LineWidth', 2);   %starts at day 1
set(p, 'Color', [0.5 0 0.5]);   %yellow doesn't show up
hold on;

%% mouse Vgateight: randomized

trainingDay3 = num(timespan3, 5);
percCorrectsTot3 = num(timespan3, 22);

percCorrectsMat3(:,1) = trainingDay3;
percCorrectsMat3(:,2) = percCorrectsTot3;

logMat3 = isfinite(percCorrectsMat3);

for i3 = 1:length(percCorrectsMat3)
    if logMat3(i3,1) == 1 && logMat3(i3,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat3(i3,:) = percCorrectsMat3(i3,:);
    end
end

vecTrain3 = find(percCorrectsTotMat3(:,1) > 0);
vecPerform3 = find(percCorrectsTotMat3(:,2) > 0);
performMat3(:,1) = percCorrectsTotMat3(vecTrain3,1);
performMat3(:,2) = percCorrectsTotMat3(vecPerform3,2);

vec3one = performMat3(:,1); vec3two = performMat3(:,2);
time2perform3one1 = find(vec3two == 75);

time2perform3one = vec3one(time2perform3one1);

plot(performMat3(:,1), performMat3(:,2), 'm', 'LineWidth', 2);   %starts at day 3
xlabel('Training day', 'FontSize', 28); xlim([0 90]);
ylabel('Performance (% corrects)','FontSize', 28); ylim([30 100]);
        ax = gca; 
        ax.FontSize = 28;
% title(mousename3);

%% mouse SC2:
trainingDay5 = num(timespan5, 5);
percCorrectsTot5 = num(timespan5, 22);

trainingDay4 = num(timespan4, 5);
percCorrectsTot4 = num(timespan4, 22);

percCorrectsMat5(:,1) = trainingDay5;
percCorrectsMat5(:,2) = percCorrectsTot5;

logMat5 = isfinite(percCorrectsMat5);
for i3 = 1:length(percCorrectsMat5)
    if logMat5(i3,1) == 1 && logMat5(i3,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat5(i3,:) = percCorrectsMat5(i3,:);
    end
end
vecTrain5 = find(percCorrectsTotMat5(:,1) > 0);
vecPerform5 = find(percCorrectsTotMat5(:,2) > 0);
performMat5(:,1) = percCorrectsTotMat5(vecTrain5,1);
performMat5(:,2) = percCorrectsTotMat5(vecPerform5,2);

vec5one = performMat5(:,1); vec5two = performMat5(:,2);
time2perform5one1 = find(vec5two == 75);

time2perform5one = vec5one(time2perform5one1);
    ax = gca;
%     mn = min(performMat5(:,1)); 
    mx = max(performMat5(:,1));
    y1 = 75; y2 = 75;
    l = line([0 mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat5(6:75,1), performMat5(6:75,2), 'r', 'LineWidth', 2);  %starts at day 4

ylabel('Performance (% corrects)', 'FontSize', 28); ylim([30 100]);
xlabel('Training day', 'FontSize', 28); xlim([0 90]);
        ax = gca; 
        ax.FontSize = 28;

%% Calculate the mean of randomized sessions' performances:
performMatNew2_2 = performMat1two(:,2); 
performMatNew2_3 = performMat1two(:,3);
performMat4avg2 = performMatNew2_2;

performMatNew3_2 = performMat3(:,2);
performMatNew3_1 = performMat3(:,1);
performMat4avg3 = performMatNew3_2;

performMatNew5_2 = performMat5(:,2);
performMatNew5_1 = performMat5(:,1);
performMat4avg5 = performMatNew5_2; % Align performance vecs from the last session:

performMat = NaN(length(performMat4avg5), 3);
for k = 1:length(performMat4avg2)
    performMat(k,1) = performMat4avg2(k);
end
for k = 1:length(performMat4avg3)
    performMat(k,2) = performMat4avg3(k);
end
for k = 1:length(performMat4avg5);
    performMat(k,3) = performMat4avg5(k);
end

% %now that matrix is filled from back-end, calculate mean on front-end:
% performMat = flip(performMat,1);
% performRandAvg = nanmean(performMat,2);
% plot(performRandAvg, 'k', 'lineWidth', 4);
% hold off;       
 
% calculate mean performance: better at finding the in-betweens compared to median
performRandAvg = nanmean(performMat,2);
% plot(performRandAvg, 'k', 'lineWidth', 4);

% shift where the mean is plotted based on the average x values of individual vector starts
x_startMeanRand = round((1+3+8)/3);
xminRand = x_startMeanRand;
xmaxRand = length(performMat4avg5) + (x_startMeanRand-1);
plot([xminRand:xmaxRand], performRandAvg, 'k-', 'lineWidth', 4);

hold off;       
       
% calculat median
% performRandMedian = nanmedian(performMat,2);
% plot(performRandMedian, 'k', 'lineWidth', 4);
% hold off;       

             
%% mouse Vgatseven - blocks

trainingDay2 = num(timespan2, 5);
percCorrectsTot2 = num(timespan2, 22);

trainingDay2two = num(timespan2two, 5);
percCorrectsTot2two = num(timespan2two, 22);

percCorrectsMat2(:,1) = trainingDay2;
percCorrectsMat2(:,2) = percCorrectsTot2;

percCorrectsMat2two(:,1) = trainingDay2two;
percCorrectsMat2two(:,2) = percCorrectsTot2two;

logMat2 = isfinite(percCorrectsMat2);
logMat2two = isfinite(percCorrectsMat2two);

for i2 = 1:length(percCorrectsMat2)
    if logMat2(i2,1) == 1 && logMat2(i2,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat2(i2,:) = percCorrectsMat2(i2,:);
    end
end
for i3 = 1:length(percCorrectsMat2two)
    if logMat2two(i3,1) == 1 && logMat2two(i3,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat2two(i3,:) = percCorrectsMat2two(i3,:);
    end
end

vecTrain2 = find(percCorrectsTotMat2(:,1) > 0);
vecPerform2 = find(percCorrectsTotMat2(:,2) > 0);
performMat2(:,1) = percCorrectsTotMat2(vecTrain2,1);
performMat2(:,2) = percCorrectsTotMat2(vecPerform2,2);

vecTrain2two = find(percCorrectsTotMat2two(:,1) > 0);
vecPerform2two = find(percCorrectsTotMat2two(:,2) > 0);
performMat2two(:,1) = percCorrectsTotMat2two(vecTrain2two,1);
performMat2two(:,2) = percCorrectsTotMat2two(vecPerform2two,2);
performMat2two(:,3) = 1:length(performMat2two(:,2));    %normalize sessions for new coupling

vec2one = performMat2(:,1); vec2two = performMat1(:,2);
vec2three = performMat2two(:,1); vec2four = performMat2two(:,2);
time2perform2one1 = find(vec2two == 75);
time2perform2two1 = find(vec2four == 75);

time2perform2one = vec2one(time2perform2one1);
time2perform2two = vec2three(time2perform2two1(1));

% Plot 2
figure; hold on;
%     ax = gca;
%     mn = min(performMat2(:,1)); mx = max(performMat2two(:,1));
%     y1 = 75; y2 = 75;
%     l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
%     set(l,'LineStyle','--')

% plot(performMat2(:,1), performMat2(:,2), 'cyan', 'LineWidth', 2);
plot(performMat2two(4:end,3), performMat2two(4:end,2), 'bl', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 28); ylim([30 100]); %xlim([8 max(performMat2two(:,1))]);
xlabel('Training day', 'FontSize', 28); xlim([0 90]);
% title(mousename2);        
ax = gca; 
ax.FontSize = 28;

%% mouse SC1:

trainingDay4 = num(timespan4, 5);
percCorrectsTot4 = num(timespan4, 22);

percCorrectsMat4(:,1) = trainingDay4;
percCorrectsMat4(:,2) = percCorrectsTot4;

logMat4 = isfinite(percCorrectsMat4);

for i3 = 1:length(percCorrectsMat4)
    if logMat4(i3,1) == 1 && logMat4(i3,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat4(i3,:) = percCorrectsMat4(i3,:);
    end
end

vecTrain4 = find(percCorrectsTotMat4(:,1) > 0);
vecPerform4 = find(percCorrectsTotMat4(:,2) > 0);
performMat4(:,1) = percCorrectsTotMat4(vecTrain4,1);
performMat4(:,2) = percCorrectsTotMat4(vecPerform4,2);

vec4one = performMat4(:,1); vec4two = performMat4(:,2);
time2perform4one1 = find(vec4two == 75);

time2perform4one = vec4one(time2perform4one1);

% Plot 4
% figure; hold on;
%     ax = gca;
%     mn = min(performMat4(:,1)); mx = max(performMat4(:,1));
%     y1 = 75; y2 = 75;
%     l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
%     set(l,'LineStyle','--')

plot(performMat4(:,1), performMat4(:,2), 'cyan', 'LineWidth', 2);

% ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]);
% xlabel('Training day', 'FontSize', 18); xlim([0 90]);
% title(mousename4);


%% mouse SC3:

trainingDay6 = num(timespan6, 5);
percCorrectsTot6 = num(timespan6, 22);

percCorrectsMat6(:,1) = trainingDay6;
percCorrectsMat6(:,2) = percCorrectsTot6;

logMat6 = isfinite(percCorrectsMat6);

for i6 = 1:length(percCorrectsMat6)
    if logMat6(i6,1) == 1 && logMat6(i6,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat6(i6,:) = percCorrectsMat6(i6,:);
    end
end

vecTrain6 = find(percCorrectsTotMat6(:,1) > 0);
vecPerform6 = find(percCorrectsTotMat6(:,2) > 0);
performMat6(:,1) = percCorrectsTotMat6(vecTrain6,1);
performMat6(:,2) = percCorrectsTotMat6(vecPerform6,2);

vec6one = performMat6(:,1); vec6two = performMat6(:,2);
time2perform6one1 = find(vec6two == 75);
time2perform6one = vec6one(time2perform6one1);

% Plot 6
% figure; hold on;
    ax = gca;
%     mn = min(performMat6(:,1)); 
    mx = max(performMat6(:,1));
    y1 = 75; y2 = 75;
    l = line([0 mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat6(:,1), performMat6(:,2), 'g', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 28); ylim([30 100]);
xlabel('Training day', 'FontSize', 28); xlim([0 30]);
% title(mousename6);
        ax = gca; 
        ax.FontSize = 28;


%% Calculate the mean of blocked sessions' performances:

performMatNew2two_2 = performMat2two(:,2); %correct for different sized vectors & diff training day starts
performMatNew2two_3 = performMat2two(:,3);
performMat4avg2two = performMatNew2two_2;

performMatNew4_1 = performMat4(:,1);
performMatNew4_2 = performMat4(:,2);
performMat4avg4 = performMatNew4_2;

performMatNew6_2 = performMat6(:,2);
performMatNew6_1 = performMat6(:,1);
performMat4avg6 = performMatNew6_2;

performMatBlk = NaN(length(performMat4avg2two), 3);
for k = 1:length(performMat4avg2two)
    performMatBlk(k,1) = performMat4avg2two(k);
end
for k = 1:length(performMat4avg4)
    performMatBlk(k,2) = performMat4avg4(k);
end
for k = 1:length(performMat4avg6)
    performMatBlk(k,3) = performMat4avg6(k);
end

performBlkAvg = nanmean(performMatBlk,2);
% plot(performBlkAvg, 'k', 'lineWidth', 2);

% shift where the mean is plotted based on the average x values of individual vector starts
x_startMeanBlk = round((4+8)/2);  %2 vecs start at 4 & 1 at 8, but I don't want 4 weighted
xminBlk = x_startMeanBlk;
xmaxBlk = length(performMatBlk) + (x_startMeanBlk-1);
plot([xminBlk:xmaxBlk], performBlkAvg, 'k-', 'lineWidth', 4);

% performBlkMedian = nanmedian(performMatBlk,2);
% plot(performBlkMedian, 'r', 'lineWidth', 1)

%% Try something else to get mean of diff sized vecs :
% vec1 = [0 9 5 1 1 1 1 1];
% vec2 = [2 8 7 0 0 0 0 0];
% vec3 = [5 4 3 2 4 6 0 0];
% 
% dataset = vertcat(vec1, vec2, vec3);
% dataset = dataset';
% 
% f = dataset(:,:,1);
%  d = diff(dataset,1,3);
% %  get rid of the negatives (assuming your data is purely monotonic
%  d(d<0) = 0; 
% %  Put them together
%  d =  f(:,:,2:(size(d,3)+1));
% %  Take the cumulative sum:
%  f = cumsum(f,3);
% %  means = mean(f,chosen_dim);
%  means = mean(f,2);
% 
%         
%% Set up vecs for plotting number of days to reach 75%

% Vgattwo = [2 150];
% Vgatthree = [2 160];
% Vgatfour =  [2 time2perform(1)];        %won't have a second vec
% Vgatfive1 = [2 92];                     %no, he stayed w/ first coupling  %since the first vec was empty, training day = 90+
% % Vgatfive2 = [8 time2perform1two(1)];
% Vgatsix = [2 89];
% Vgatseven1 = [2 90];                     %since the first vec was empty, training day = 90+
% % Vgatseven2 =[4 time2perform2two(1)];
% Vgatseven2 = [4 8]; %since 8 d after switch, he was > 75%; 
% Vgatseven = vertcat(Vgatseven1, Vgatseven2);
% 
% Vgateight2 = [4 time2perform3one(1)];
% % SC1 = [2 time2perform4one];
% SC1 = [4 12];
% % SC2 = [4 time2perform5one(1)];
% SC2 = [4 20];
% 
% %Avg. the blocks for old regime and compare to new; same for randomized:
% blocksVec1 = vertcat(Vgattwo(2), Vgatfour(2), Vgatsix (2), Vgatseven1(2));
% blocksMean1 = mean(blocksVec1); blocksstd1 = std(blocksVec1); semBlocks1 = blocksstd1/sqrt(4);
% 
% blocksVec2 = vertcat(Vgatseven2(2), SC1(2));
% blocksMean2 = mean(blocksVec2); blocksstd2 = std(blocksVec2); semBlocks2 = blocksstd1/sqrt(2);
% 
% randVec1 = vertcat(Vgatthree(2), Vgatfive1(2));
% randMean1 = mean(randVec1); randstd1 = std(randVec1); semRand1 = randstd1/sqrt(2);
% 
% randVec2 = vertcat(Vgateight2(2), SC2(2));
% randMean2 = mean(randVec2); randstd2 = std(randVec2); semRand2 = randstd2/sqrt(2);
% 
% %% final plotting:
% 
% figure; hold on;
% ylabel('# Training days to 75% correct', 'FontSize', 18);    %ylim([0 100]);
% ylim([0 140]); xlim([0 8]);
% 
% h = plot(2,blocksMean1, 'bl.', 'MarkerSize', 24);
% i = plot(4,blocksMean2, 'bl.', 'MarkerSize', 24);
% j = plot(2.2,randMean1, 'r.', 'MarkerSize', 24);
% k = plot(4.2,randMean2, 'r.', 'MarkerSize', 24);
% 
% l = errorbar(2,blocksMean1,semBlocks1, 'bl-', 'LineWidth', 2);
% m = errorbar(4,blocksMean2,semBlocks2, 'bl-', 'LineWidth', 2);
% n = errorbar(2.2,randMean1,semRand1, 'r-', 'LineWidth', 2);
% o = errorbar(4.2,randMean2,semRand2, 'r-', 'LineWidth', 2);
% 
% set(gca, 'YTick', 0:20:140, 'XTick', 0:1:4);
% ax = gca; 
% set(ax);
% ax.XTickLabelMode = 'manual';
% ax.XTickLabel = {'', '','Old', '', 'New'};
% ax.FontSize = 18;
% 
