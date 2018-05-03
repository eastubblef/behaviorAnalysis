function [theCatchTrials] = catchTrials(taskbase)

% Purpose of this function is to extract catch trials from block sessions
% 7.12.17
% Also pulls out stim trials & determine which were catch trials
% Must first run behaveLoad6working to obtain the taskbase structure file (*tb.mat) - this is the input file

fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC3';
filenamestr = 'SC3';

%% Load taskbase structure for behavioral data 
                      
if exist('filenamestr', 'var');
    [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       
end
taskbase = strcat(path, filenamestrT);
load(taskbase);

trialMat = taskbase.trialMat;
trialMat14 = trialMat(:,14);
trialMat1 = trialMat(:,1);

%% Make all R-move trials -1 and L-move trials +1; find derivative (-2, +2); plot out for inflection-visualization

trialMat15 = [];
for i = 1:length(trialMat14)
    if trialMat14(i) < 0
        trialMat15(i) = -1;
    else if trialMat14(i) > 0
            trialMat15(i) = 1;
        end
    end
end
trialMat15 = trialMat15';
plot(trialMat15); ylim([-3 3]);
diffVec = diff(trialMat15);   %-2,+2 = catch trial; one alone = block switch

diffVecAdjust = [];
for j = 1:length(diffVec)
    diffVecAdjust(j+1) = diffVec(j);
end
diffVecAdjust = diffVecAdjust';
catchTrialPre = find(abs(diffVecAdjust) > 0);        

catchTrial1 = diff(catchTrialPre);
for k = 1:length(catchTrial1)
    catchTrialAdjust(k+1) = catchTrial1(k);
end
catchTrialAdjust = catchTrialAdjust';  %get rid of the ones 

for l = 1:length(catchTrialAdjust)
    if catchTrialAdjust(l) == 1
        catchTrialPre(l) = 0;
    end
end

catchTrialsAlmost = find(catchTrialPre > 0);
theCatchTrials = catchTrialPre(catchTrialsAlmost);


end

