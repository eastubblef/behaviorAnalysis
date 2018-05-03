
%% This script is for generating a heat map (imagesc) of binned rx times' %corrects 
% Bin sizes are set for 200 ms = fastest bin; 200-600 ms = fast bin; 600-1000 = mid bin; >1000 = slow bin
% 11.2017  This script calls rxSort for forming the matrix
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
RXmat1 = rxMat1';
clear taskbaseFile1;

load(taskbaseFile2);
rxMat2 = taskbase.rxMat;
RXmat2 = rxMat2';
clear taskbaseFile2;

load(taskbaseFile3);
rxMat3 = taskbase.rxMat;
RXmat3 = rxMat3';
clear taskbaseFile3;

load(taskbaseFile4);
rxMat4 = taskbase.rxMat;
RXmat4 = rxMat4';
clear taskbaseFile4;

% rxVec1 = rxMat1(:,2)
% RXmat = [];
% for i = 1:length(rxVec1)
%     RXmat(:,i) = rxVec1(i);
% end

[rxSort] = rxSort4Mat(rxMat1);
rxSort1 = rxSort';
clear rxSort;

[rxSort] = rxSort4Mat(rxMat2);
rxSort2 = rxSort';
clear rxSort;

[rxSort] = rxSort4Mat(rxMat3);
rxSort3 = rxSort';
clear rxSort;

[rxSort] = rxSort4Mat(rxMat4);
rxSort4 = rxSort';
clear rxSort;

[NumFastest, NumFast, NumMid, NumSlow] = rxSortMat4perform(rxSort1);
NumFastest1 = NumFastest;
NumFast1 = NumFast;
NumMid1 = NumMid;
NumSlow1 = NumSlow;
cat1 = horzcat(NumFastest1, NumFast1, NumMid1, NumSlow1);

[NumFastest, NumFast, NumMid, NumSlow] = rxSortMat4perform(rxSort2);
NumFastest2 = NumFastest;
NumFast2 = NumFast;
NumMid2 = NumMid;
NumSlow2 = NumSlow;
cat2 = horzcat(NumFastest2, NumFast2, NumMid2, NumSlow2);

[NumFastest, NumFast, NumMid, NumSlow] = rxSortMat4perform(rxSort3);
NumFastest3 = NumFastest;
NumFast3 = NumFast;
NumMid3 = NumMid;
NumSlow3 = NumSlow;
cat3 = horzcat(NumFastest3, NumFast3, NumMid3, NumSlow3);

[NumFastest, NumFast, NumMid, NumSlow] = rxSortMat4perform(rxSort4);
NumFastest4 = NumFastest;
NumFast4 = NumFast;
NumMid4 = NumMid;
NumSlow4 = NumSlow;

%% Create bins for rx times and insert % corrects
catPerform = vertcat(cat1,cat2,cat3,cat4);

h_image = imagesc(catPerform);



