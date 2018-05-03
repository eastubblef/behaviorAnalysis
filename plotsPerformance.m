%% This script plots total percent performance for a behavioral session.
%Just save your excel file as excel 5.0/95 then it will work - Mathworks
% PLOTS EA. MOUSE SEPARATELY
%blocks + full bar (old coupling = solid) = cyan; blocks + half bar (new coupling = dotted) = dk blue; 
%randomized + full bar (old coupling = solid) = magenta; randomized + half bar (new coupling = dotted) = red

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
timespan3 = 981:1103                  %update; includes blocks, starting out, then conversion to randomized
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

% %mouse Vgatfour: OLD COUPLING
% trainingDay = num(timespan, 1);
% percCorrectsTot = num(timespan, 18);
% 
% percCorrectsMat(:,1) = trainingDay;
% percCorrectsMat(:,2) = percCorrectsTot;
% 
% logMat = isfinite(percCorrectsMat);
% 
% for i = 1:length(percCorrectsMat)
%     if logMat(i,1) == 1 && logMat(i,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
%         percCorrectsTotMat(i,:) = percCorrectsMat(i,:);
%     end
% end
% 
% vecTrain = find(percCorrectsTotMat(:,1) > 0);
% vecPerform = find(percCorrectsTotMat(:,2) > 0);
% performMat(:,1) = percCorrectsTotMat(vecTrain,1)
% performMat(:,2) = percCorrectsTotMat(vecPerform,2);
% 
% vec1 = performMat(:,1); vec2 = performMat(:,2);
% time2perform1 = find(vec2 == 75);
% time2perform = vec1(time2perform1);
% 
% % Plot like Fig 1 of Burgess et al., 2016 
% figure; hold on;
%     ax = gca;
% %     mn = min(performMat(:,1)); mx = max(performMat(:,1));
%     mn = 10; mx = max(performMat(:,1));
%     y1 = 75; y2 = 75;
%     l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
%     set(l,'LineStyle','--')
% 
% plot(performMat(:,1), performMat(:,2), 'c', 'LineWidth', 2);
% ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]);
% xlabel('Training day', 'FontSize', 18); 
% title(mousename);
% hold off;

%% mouse Vgatfive

% trainingDay1 = num(timespan1, 1);
% percCorrectsTot1 = num(timespan1, 18);

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
performMat1(:,1) = percCorrectsTotMat1(vecTrain1,1)
performMat1(:,2) = percCorrectsTotMat1(vecPerform1,2);

vecTrain1two = find(percCorrectsTotMat1two(:,1) > 0);
vecPerform1two = find(percCorrectsTotMat1two(:,2) > 0);
performMat1two(:,1) = percCorrectsTotMat1two(vecTrain1two,1)
performMat1two(:,2) = percCorrectsTotMat1two(vecPerform1two,2);

vec1one = performMat1(:,1); vec1two = performMat1(:,2);
vec1three = performMat1two(:,1); vec1four = performMat1two(:,2);
time2perform1one1 = find(vec1two == 75);
time2perform1two1 = find(vec1four == 75);

time2perform1one = vec1one(time2perform1one1);
time2perform1two = vec1three(time2perform1two1);

% Plot 
figure; hold on;
    ax = gca;
    mn = min(performMat1(:,1)); mx = max(performMat1two(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat1(:,1), performMat1(:,2), 'm', 'LineWidth', 2);
plot(performMat1two(:,1), performMat1two(:,2), 'r', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]); xlim([min(performMat1(:,1)) max(performMat1two(:,1))]);
xlabel('Training day', 'FontSize', 18); 
title(mousename1);

%% mouse Vgatseven

% trainingDay2 = num(timespan2, 1);
% percCorrectsTot2 = num(timespan2, 18);
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
performMat2(:,1) = percCorrectsTotMat2(vecTrain2,1)
performMat2(:,2) = percCorrectsTotMat2(vecPerform2,2);

vecTrain2two = find(percCorrectsTotMat2two(:,1) > 0);
vecPerform2two = find(percCorrectsTotMat2two(:,2) > 0);
performMat2two(:,1) = percCorrectsTotMat2two(vecTrain2two,1)
performMat2two(:,2) = percCorrectsTotMat2two(vecPerform2two,2);

vec2one = performMat2(:,1); vec2two = performMat1(:,2);
vec2three = performMat2two(:,1); vec2four = performMat2two(:,2);
time2perform2one1 = find(vec2two == 75);
time2perform2two1 = find(vec2four == 75);

time2perform2one = vec2one(time2perform2one1);
time2perform2two = vec2three(time2perform2two1(1));

% Plot 2
figure; hold on;
    ax = gca;
    mn = min(performMat2(:,1)); mx = max(performMat2two(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat2(:,1), performMat2(:,2), 'cyan', 'LineWidth', 2);
plot(performMat2two(:,1), performMat2two(:,2), 'bl', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]); xlim([8 max(performMat2two(:,1))]);
xlabel('Training day', 'FontSize', 18); xlim([0 90]);
title(mousename2);

%% mouse Vgateight:

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
performMat3(:,1) = percCorrectsTotMat3(vecTrain3,1)
performMat3(:,2) = percCorrectsTotMat3(vecPerform3,2);

vec3one = performMat3(:,1); vec3two = performMat3(:,2);
time2perform3one1 = find(vec3two == 75);

time2perform3one = vec3one(time2perform3one1);

% Plot 3
figure; hold on;
    ax = gca;
    mn = min(performMat3(:,1)); mx = max(performMat3(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat3(:,1), performMat3(:,2), 'r', 'LineWidth', 2);

ylabel('Performance (% corrects)','FontSize', 18); ylim([30 100]); xlim([8 max(performMat3(:,1))]);
xlabel('Training day', 'FontSize', 18); xlim([0 90]);
title(mousename3);

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
performMat4(:,1) = percCorrectsTotMat4(vecTrain4,1)
performMat4(:,2) = percCorrectsTotMat4(vecPerform4,2);

vec4one = performMat4(:,1); vec4two = performMat4(:,2);
time2perform4one1 = find(vec4two == 75);

time2perform4one = vec4one(time2perform4one1);

% Plot 4
figure; hold on;
    ax = gca;
    mn = min(performMat4(:,1)); mx = max(performMat4(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat4(:,1), performMat4(:,2), 'bl', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]);
xlabel('Training day', 'FontSize', 18); xlim([0 90]);
title(mousename4);

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
performMat5(:,1) = percCorrectsTotMat5(vecTrain5,1)
performMat5(:,2) = percCorrectsTotMat5(vecPerform5,2);

vec5one = performMat5(:,1); vec5two = performMat5(:,2);
time2perform5one1 = find(vec5two == 75);

time2perform5one = vec5one(time2perform5one1);

% Plot 5
figure; hold on;
    ax = gca;
    mn = min(performMat5(:,1)); mx = max(performMat5(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

% plot(performMat5(5:8,1), performMat5(5:8,2), 'r', 'LineWidth', 2);
% plot(performMat5(1:5,1), performMat5(1:5,2), 'bl', 'LineWidth', 2);
plot(performMat5(6:75,1), performMat5(6:75,2), 'r', 'LineWidth', 2);
% plot(performMat5(1:5,1), performMat5(1:5,2), 'bl', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]);
xlabel('Training day', 'FontSize', 18); xlim([0 90]);

title(mousename5);

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
performMat6(:,1) = percCorrectsTotMat6(vecTrain6,1)
performMat6(:,2) = percCorrectsTotMat6(vecPerform6,2);

vec6one = performMat6(:,1); vec6two = performMat6(:,2);
time2perform6one1 = find(vec6two == 75);
time2perform6one = vec6one(time2perform6one1);

% Plot 6
figure; hold on;
    ax = gca;
    mn = min(performMat6(:,1)); mx = max(performMat6(:,1));
    y1 = 75; y2 = 75;
    l = line([mn mx], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

plot(performMat6(:,1), performMat6(:,2), 'bl', 'LineWidth', 2);

ylabel('Performance (% corrects)', 'FontSize', 18); ylim([30 100]);
xlabel('Training day', 'FontSize', 18); xlim([0 30]);
title(mousename6);


%% Set up vecs for plotting number of days to reach 75%

Vgattwo = [2 150];
Vgatthree = [2 160];
Vgatfour =  [2 time2perform(1)];        %won't have a second vec
Vgatfive1 = [2 92];                     %no, he stayed w/ first coupling  %since the first vec was empty, training day = 90+
% Vgatfive2 = [8 time2perform1two(1)];
Vgatsix = [2 89];
Vgatseven1 = [2 90];                     %since the first vec was empty, training day = 90+
% Vgatseven2 =[4 time2perform2two(1)];
Vgatseven2 = [4 8]; %since 8 d after switch, he was > 75%; 
Vgatseven = vertcat(Vgatseven1, Vgatseven2);

Vgateight2 = [4 time2perform3one(1)];
% SC1 = [2 time2perform4one];
SC1 = [4 12];
% SC2 = [4 time2perform5one(1)];
SC2 = [4 20];

%Avg. the blocks for old regime and compare to new; same for randomized:
blocksVec1 = vertcat(Vgattwo(2), Vgatfour(2), Vgatsix (2), Vgatseven1(2));
blocksMean1 = mean(blocksVec1); blocksstd1 = std(blocksVec1); semBlocks1 = blocksstd1/sqrt(4);

blocksVec2 = vertcat(Vgatseven2(2), SC1(2));
blocksMean2 = mean(blocksVec2); blocksstd2 = std(blocksVec2); semBlocks2 = blocksstd1/sqrt(2);

randVec1 = vertcat(Vgatthree(2), Vgatfive1(2));
randMean1 = mean(randVec1); randstd1 = std(randVec1); semRand1 = randstd1/sqrt(2);

randVec2 = vertcat(Vgateight2(2), SC2(2));
randMean2 = mean(randVec2); randstd2 = std(randVec2); semRand2 = randstd2/sqrt(2);

%% final plotting:

figure; hold on;
ylabel('# Training days to 75% correct', 'FontSize', 18);    %ylim([0 100]);
ylim([0 140]); xlim([0 8]);

h = plot(2,blocksMean1, 'bl.', 'MarkerSize', 24);
i = plot(4,blocksMean2, 'bl.', 'MarkerSize', 24);
j = plot(2.2,randMean1, 'r.', 'MarkerSize', 24);
k = plot(4.2,randMean2, 'r.', 'MarkerSize', 24);

l = errorbar(2,blocksMean1,semBlocks1, 'bl-', 'LineWidth', 2);
m = errorbar(4,blocksMean2,semBlocks2, 'bl-', 'LineWidth', 2);
n = errorbar(2.2,randMean1,semRand1, 'r-', 'LineWidth', 2);
o = errorbar(4.2,randMean2,semRand2, 'r-', 'LineWidth', 2);

set(gca, 'YTick', 0:20:140, 'XTick', 0:1:4);
ax = gca; 
set(ax);
ax.XTickLabelMode = 'manual';
ax.XTickLabel = {'', '','Old', '', 'New'};
ax.FontSize = 18;

