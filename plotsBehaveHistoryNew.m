%% plotsBehaveHistory
% This calls the blocks and rand mfiles for calculating dependency of choice on current trial on rewarded trial history
% 
% Run this after behaveLoad5historyNew.m
% EAS 8.2016
%% Initialize 

% initPath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/151105/2ndpass/better/behaveChunks';
initPath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/Vgatfive/forAnalysis';


% taskbaseFile = 'mVgatthreetag_2015_11_05_14_tb.mat';
taskbaseFile = 'mVgatfiveR_2016_08_25_15_tb.mat';

taskbase = strcat(initPath,'/', taskbaseFile);
load(taskbase);

LRtrialStruct = taskbase.LRtrialStruct;

%0th trial; current decision, or vector B
choiceVec = taskbase.choiceVec;

rewTrial = [];
for i = 1:length(LRtrialStruct)
    if LRtrialStruct(i) == choiceVec(i)
        rewTrial(i) = 1;
    else
        rewTrial(i) = 0;
    end
end

%matrix A, given the 0 1 rewarded trial structure
rewMat(1,:) = rewTrial;
rewMat(2,:) = choiceVec;

%new matrixA, given the -1 1 rewarded trial structure
rewMat2(1,:) = LRtrialStruct;
rewMat2(2,:) = choiceVec;

%test matrixA
% rewMat3(1,:) = [0 0 0 0 0 0 1 0 0 0];
% rewMat3(2,:) = [1 0 0 0 0 0 1 1 0 0]; %gives [-1e-16 1]

% rewMat3(1,:) = [1 1 1 1 1 1 -1 1 1 1]; rewarded or not
                %[rew         not     ] 
% rewMat3(2,:) = [-1 1 1 1 1 1 -1 -1 1 1]; choice (-1 is L)
                %[-2 2 2 2 2 2 -1 -2 2 2]; choice (-2 L rew; -1 L not; 2 is R rewarded
% choiceVec3 = rewMat3(2,:);  %gives [3.12e-16 1]

rewMat3(1,:) = [1 1 1 1 1 1 -1 1 1 1];
rewMat3(2,:) = [-2 2 2 2 2 2 -1 -2 2 2]; %gives [-5e-16 1]
choiceVec3 = rewMat3(2,:);


% %A*x = B
% x = pinv(rewMat);

%Try B/A = x:
newX = choiceVec/rewMat;

%Try pinv(A)B = x in terms of the 0 1 rewarded trial structure:
solveX = (pinv(rewMat))'*choiceVec'; 

%Try pinv(A)B = x in terms of the -1 v 1 rewarded trial structure:
solveX2 = (pinv(rewMat2))'*choiceVec';
solveX3 = (pinv(rewMat3))'*choiceVec3';     

%% Blocks: Plot Prob(left trial i) as a function of trial history:
figure; hold on;
% Prob_Li_jitter = Prob_Li + 0.01;
% 
% y = [Prob_Li_jitter, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1];        
% h = plot(y, 'r.', 'MarkerSize', 25, 'Color', [0.5 0.5 0.5]);
    ylabel('P choice (Left_i)', 'FontSize', 22);  
    xlabel('Trial history (i-num trials back)', 'FontSize', 22); 
    ylim([0 1.05]);
    set(gca, 'YTick', 0:.20:1.2, 'XTick', 0:1:8);
    ax = gca; 
    set(ax);
    ax.XTickLabelMode = 'manual'
    ax.XTickLabel = {' ', 'i','i-1','i-2','i-3', 'i-4','i-5','i-6','i-7'}
    ax.FontSize = 22;


%% Initialize for the randomized file:
% 
% fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarApr16/Vgatfive/forAnalysis';
% 
% %% Rand: Plot Prob(left trial i) as a function of trial history: 
% 
% gcf; 
% hold on; 
% y = [Prob_Li_jitter, Prob_LiMinusOne1, Prob_LiMinusTwo1, Prob_LiMinusThree1, Prob_LiMinusFour1, Prob_LiMinusFive1, Prob_LiMinusSix1, Prob_LiMinusSeven1];        
% h = plot(y, 'r.', 'MarkerSize', 25, 'Color', [0.5 0.5 0.5]);
% %    ylabel('P choice (Left)', 'FontSize', 20);  %Divided ea. by Prob_Li, so yes, P(L choice) on trial i
%     ylabel('P choice (Left_i)', 'FontSize', 22);  
%     xlabel('Trial history (i-num trials back)', 'FontSize', 22); 
%     ylim([0 1.05]);
%     set(gca, 'YTick', 0:.20:1.2, 'XTick', 0:1:8);
%     ax = gca; 
%     set(ax);
%     ax.XTickLabelMode = 'manual'
%     ax.XTickLabel = {' ', 'i','i-1','i-2','i-3', 'i-4','i-5','i-6','i-7'}
%     ax.FontSize = 22;
% 
% 
