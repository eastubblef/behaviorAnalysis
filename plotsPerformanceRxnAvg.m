%% This script plots rxn time average across behavioral sessions.
%Just save your excel file as excel 5.0/95 then it will work - Mathworks
%PLOTS ALL MICE TOGETHER for new coupling only!

%Updated 1.17.18 for better mean overlay 

fileDirectory = '/Users/stubblefielde/Desktop/mfiles/behavior/';
% cd fileDirectory;
% file = 'mouseBehavior2.xls';        %updated 9.15.17
file = 'mouseBehavior4.xls';        %updated 3.30.18
% file = 'mouseBehavior5.xls';        %updated 4.2.18

filename = strcat(fileDirectory, file);

mousename = 'Vgatfour';             %input the relevant rows from xls spreadsheet
timespan = 256:329;                 %updated 9.15.17
blocks = 1;
fullbar = 1;
newCoupling = 0;                    %drk blue - old coupling

%% Begin new coupling:
mousename1 = 'Vgatfive';
% timespan1two = 529:591;             %half bar up to recording day2
% timespan1two = 428:596;             %half bar up to recording day2; works 3.6.18
timespan1two = 345:596;                 %3.26.18
blocks1two = 0;
fullbar1two = 0;
newCoupling1two = 0; %purple

mousename2 = 'Vgatseven'; 
% timespan2 = 688:731                 %full bar, old coupling
% timespan2 = 800:856                 %full bar, old coupling
% blocks2 = 1;
% fullbar2 = 1,                       %solid
% timespan2two = 742:824;             %half bar, new coupling
timespan2two = 820:965;
blocks2two = 1;
fullbar2two = 0;
newCoupling2two = 1; %cyan

mousename3 = 'Vgateight';
% timespan3 = 838:911                 %includes blocks, starting out, then conversion to randomized
timespan3 = 981:1104;                  %update; includes blocks, starting out, then conversion to randomized
blocks3 = 0;
fullbar3 = 0;
newCoupling3 = 1;   %magenta

mousename4 = 'SC1';
% timespan4 = 926:940;
% timespan4 = 1116:1222;                %update
timespan4 = 1116:1266;                %update
blocks4 = 1;
fullbar4 = 0;
newCoupling4 = 1;   %cyan

mousename5 = 'SC2';
% timespan5 = 953:963;                   %953:967 for full blocks timespan
% timespan5 = 1278:1425;                   %update
blocks5 = 0;
fullbar5 = 0;
newCoupling5 = 1;   %cyan

timespan5two = 1282:1425;                   %update
blocks5two = 0;
fullbar5two = 0;
newCoupling5two = 1;  %magenta

mousename6 = 'SC3';
timespan6 = 1457:1503;                   %update
timespan6 = 1466:1503;                   %update

blocks6 = 1;
fullbar6 = 0;
newCoupling6 = 1;   %cyan

%% Make it an xlsRead function: 

[num, txt] = xlsread(filename); 

%% Plot like Fig 1 of Burgess et al., 2016 
% mouse Vgatfive randomized - given his rxn times, I may not be able to use
% him.

trainingDay1two = num(timespan1two, 19);
% percCorrectsTot1two = num(timespan1two, 18);
rxMat1two = num(timespan1two, 20);                      %updated for the rxn v. training day plot
percCorrectsMat1two(:,1) = trainingDay1two;
percCorrectsMat1two(:,2) = rxMat1two;

logMat1two = isfinite(percCorrectsMat1two);
for i2 = 1:length(percCorrectsMat1two)
    if logMat1two(i2,1) == 1 && logMat1two(i2,2) == 1   %eliminate any row in which a Nan appears to keep consistent between cols
        percCorrectsTotMat1two(i2,:) = percCorrectsMat1two(i2,:);
    end
end

vecTrain1two = find(percCorrectsTotMat1two(:,1) > 0);
vecPerform1two = find(percCorrectsTotMat1two(:,2) > 0);
performMat1two(:,1) = percCorrectsTotMat1two(vecTrain1two,1);
performMat1two(:,2) = percCorrectsTotMat1two(vecPerform1two,2);
performMat1two(:,3) = 1:length(performMat1two(:,2));

figure; hold on;  %Vgatfive
% p = plot(performMat1two(:,3), performMat1two(:,2), 'y', 'LineWidth', 2);   %starts at day 1
p = plot(performMat1two(:,3), performMat1two(:,2), 'y.', 'MarkerSize', 24);   %starts at day 1

set(p, 'Color', [0.5 0 0.5]);   %yellow doesn't show up
hold on;

%% mouse Vgateight: randomized

trainingDay3 = num(timespan3, 19);
% percCorrectsTot3 = num(timespan3, 18);
rxMat3 = num(timespan3, 20);                      %updated for the rxn v. training day plot
percCorrectsTot3 = rxMat3;

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

% vec3one = performMat3(:,1); vec3two = performMat3(:,2);
% time2perform3one1 = find(vec3two == 75);
% time2perform3one = vec3one(time2perform3one1);

% plot(performMat3(:,1), performMat3(:,2), 'm', 'LineWidth', 2);   %starts at day 3
plot(performMat3(:,1), performMat3(:,2), 'm.', 'MarkerSize', 24);   %starts at day 3

xlabel('Training day', 'FontSize', 28); xlim([0 90]);
ylabel('Reaction times (ms)','FontSize', 28); %ylim([30 100]);
        ax = gca; 
        ax.FontSize = 28;
% title(mousename3);

%% mouse SC2:
trainingDay5 = num(timespan5two, 19);
% percCorrectsTot5 = num(timespan5two, 18);
% trainingDay4 = num(timespan4, 5);
% percCorrectsTot4 = num(timespan4, 22);
rxn5 = num(timespan5two, 20);
percCorrectsTot5 = rxn5;

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

% vec5one = performMat5(:,1); vec5two = performMat5(:,2);
% time2perform5one1 = find(vec5two == 75);
% time2perform5one = vec5one(time2perform5one1);
    ax = gca;
    mx = max(performMat5(:,1));
% plot(performMat5(:,1), performMat5(:,2), 'r', 'LineWidth', 2);  
plot(performMat5(:,1), performMat5(:,2), 'r.', 'MarkerSize', 24);  

ylabel('Reaction time (ms)', 'FontSize', 28); %ylim([30 100]);
xlabel('Training day', 'FontSize', 28); xlim([0 80]);
        ax = gca; 
        ax.FontSize = 28;

%% Calculate the mean of randomized sessions' performances:
performMatNew2_2 = performMat1two(:,2);   %Vgatfive
performMatNew2_3 = performMat1two(:,3);
performMat4avg2 = performMatNew2_2;

performMatNew3_2 = performMat3(:,2);      %Vgateight
performMatNew3_1 = performMat3(:,1);
performMat4avg3 = performMatNew3_2;

performMatNew5_2 = performMat5(:,2);      %SC2
performMatNew5_1 = performMat5(:,1);
performMat4avg5 = performMatNew5_2;       %Align performance vecs from the last session:

performMat = NaN(length(performMat4avg2), 3);
% performMat = NaN(length(performMatNew5_2), 3);

for k = 1:length(performMat4avg2)
    performMat(k,1) = NaN;
end
for k = 1:length(performMat4avg2)          %Include for Vgatfive's reaction times
    performMat(k,1) = performMat4avg2(k);
end

for k = 1:length(performMat4avg3)
    performMat(k,2) = performMat4avg3(k);
end
for k = 1:length(performMat4avg5)
    performMat(k,3) = performMat4avg5(k);
end

numelMat4size = length(performMat4avg2);
% numelMat4size = length(performMatNew5_2);  %whichever is longest
performMat = performMat(1:numelMat4size, :);

% calculate mean performance: better at finding the in-betweens compared to median
performRandAvg = nanmean(performMat,2);
plot(performRandAvg, 'k', 'lineWidth', 4);
plot(performRandAvg, 'k.', 'MarkerSize', 24);

ylim([0 6000]);
hold off;       
                   
%% mouse Vgatseven - blocks

trainingDay2two = num(timespan2two, 19);
% percCorrectsTot2two = num(timespan2two, 18);
rxnTot2two = num(timespan2two, 20);
percCorrectsTot2two = rxnTot2two;

percCorrectsMat2(:,1) = trainingDay2two;
percCorrectsMat2(:,2) = percCorrectsTot2two;

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

vecTrain2two = find(percCorrectsTotMat2two(:,1) > 0);
vecPerform2two = find(percCorrectsTotMat2two(:,2) > 0);
performMat2two(:,1) = percCorrectsTotMat2two(vecTrain2two,1);
performMat2two(:,2) = percCorrectsTotMat2two(vecPerform2two,2);
performMat2two(:,3) = 1:length(performMat2two(:,2));    %normalize sessions for new coupling

vec2four = performMat2two(:,2);
% Plot 2
figure; hold on;
plot(performMat2two(:,3), performMat2two(:,2), 'bl', 'LineWidth', 2); %Vgatseven

ylabel('Reaction time (ms)', 'FontSize', 28); ylim([0 6000]);
xlabel('Training day', 'FontSize', 28); xlim([0 80]);
% title(mousename2);        
ax = gca; 
ax.FontSize = 28;

%% mouse SC1:

trainingDay4 = num(timespan4, 19);
% percCorrectsTot4 = num(timespan4, 18);
rxnTot4 = num(timespan4, 20);
percCorrectsTot4 = rxnTot4;

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
plot(performMat4(:,1), performMat4(:,2), 'cyan', 'LineWidth', 2);   %SC1
% title(mousename4);


%% mouse SC3:

trainingDay6 = num(timespan6, 19);
rxnTot6 = num(timespan6, 20);
percCorrectsTot6 = rxnTot6;

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

% Plot 6
ax = gca;
plot(performMat6(:,1), performMat6(:,2), 'g', 'LineWidth', 2);   %SC3
ylabel('Reaction time (ms)', 'FontSize', 28); ylim([0 6000]);
xlabel('Training day', 'FontSize', 28); xlim([0 80]);
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


performMatBlk = NaN(length(performMatNew2_2), 3);
% performMat = NaN(length(performMatNew5_2), 3);

for k = 1:length(performMatNew2_2)
    performMatBlk(k,1) = NaN;
end


for k = 1:length(performMat4avg2two)
    performMatBlk(k,1) = performMat4avg2two(k);
end
for k = 1:length(performMat4avg4)
    performMatBlk(k,2) = performMat4avg4(k);
end
for k = 1:length(performMat4avg6)
    performMatBlk(k,3) = performMat4avg6(k);
end

numelMat4sizeBlk = length(performMatNew2_2);
performMatBlk = performMatBlk(1:numelMat4sizeBlk, :);

performBlkAvg = nanmean(performMatBlk,2);
plot(performBlkAvg, 'k', 'lineWidth', 4);
plot(performBlkAvg, 'k.', 'MarkerSize', 22);


%% Generate separate plots for the average rxn times across training days
%Randomized data:
figure; hold on;            %Plot randomized averages
% performRandAvg = nanmean(performMat,2);
% plot(performRandAvg, 'k', 'lineWidth', 4);
plot(performRandAvg, 'r.', 'markerSize', 24);

ylabel('Mean reaction time (ms)', 'FontSize', 28); ylim([0 5000]);
xlabel('Training day', 'FontSize', 28); xlim([0 60]);
ax = gca; 
ax.FontSize = 28;

xlim([0 60]);
ylim([0 5000]); 
hold off;   

figure; hold on;            %Plot the blocks' avg
% plot(performBlkAvg, 'k', 'lineWidth', 4);
plot(performBlkAvg, 'bl.', 'markerSize', 24);

xlabel('Training day', 'FontSize', 28); 
ylabel('Mean reaction time (ms)', 'FontSize', 28); 

ax = gca; 
ax.FontSize = 28;

xlim([0 60]);
ylim([0 5000]);
 
hold off;   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
%% Plot the average rxn times as a line:
figure; hold on;            %Plot randomized averages
% performRandAvg = nanmean(performMat,2);
plot(performRandAvg, 'r', 'lineWidth', 4);
% plot(performRandAvg, 'r.', 'markerSize', 24);

ylabel('Mean reaction time (ms)', 'FontSize', 28); ylim([0 5000]);
xlabel('Training day', 'FontSize', 28); xlim([0 60]);
ax = gca; 
ax.FontSize = 28;

xlim([0 60]);
ylim([0 5000]); 
% hold off;   
% 
% figure; hold on;            %Plot the blocks' avg
plot(performBlkAvg, 'bl', 'lineWidth', 4);
% plot(performBlkAvg, 'bl.', 'markerSize', 24);

xlabel('Training day', 'FontSize', 28); 
ylabel('Mean reaction time (ms)', 'FontSize', 28); 

ax = gca; 
ax.FontSize = 28;

xlim([0 60]);
ylim([0 5000]);
    x1 = 0; x2 = 60;
    y1 = 100; y2 = 100;
    l = line([x1 x2], [y1 y2], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

    y3 = 400; y4 = 400;
    l = line([x1 x2], [y3 y4], 'Color', 'k'); % plot line first so that data points occlude it
    set(l,'LineStyle','--')

 
hold off;   





