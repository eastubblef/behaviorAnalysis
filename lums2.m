
function [alignTimes, positions] = lums2(taskbase, contra, ipsi, lumTrials)

%Updated 8/9/17 for plotting trajectories 
%Purpose: plot positions of bar over trial duration

%Run this after behaveLoad6workingOppStim_NewBin_LvR 
%Hard code inputs: should-mv contra v. ipsi trials on line 35-36

%% Change hard-coded lumVals here:
%Pull out relavent velocities, positions, and rxTimes:
fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/SC2';
filenamestr = 'SC2';


%% Load taskbase structure for behavioral data 
                      
if exist('filenamestr', 'var');
    [filenamestrT, path] = uigetfile('*_tb.mat*','select the taskbase file', fpath)       
end
taskbase = strcat(path, filenamestrT);
load(taskbase);
mouseName = filenamestr;
numSession = 170830;

% Input luminance values - easier to find uniques, but set up the psych func's sensitivity parameters
lumVal_L1 = 255;   %should-mv left trials; ipsi for SC2
lumVal_L2 = 64;
lumVal_L3 = 15;

lumVal_R1 = -255;  %should-mv right trials; contra for SC2
lumVal_R2 = -64;
lumVal_R3 = -15;

trialMat = taskbase.trialMat;
barpositionCell = taskbase.barpositionCell;
timeTrial = taskbase.time_trial;

trialMat1 = trialMat(:,1);
trialMat4 = trialMat(:,4);
trialMat6 = trialMat(:,6);
trialMatStim = trialMat(:,13);
trialMatStartPos = trialMat(:,14);
lumTrials = trialMat(:,15);

dirLums = [];
for i = 1:length(lumTrials)
    if trialMatStartPos(i) < 0
        dirLums(i) = lumTrials(i) * -1;
    else dirLums (i) = lumTrials(i) * 1;
    end
end
trialMat(:,16) = dirLums;

% %% separate out should-mv ipsi v contra trials
% contra = find(trialMatStartPos < 1);        %for 170718; L implant - must move R to center bar
% ipsi = find(trialMatStartPos > 1);          %for 170718; L implant - must move L to center bar
LTrials = find(trialMatStartPos < 1);
RTrials = find(trialMatStartPos > 1);

% contraStimAll = trialMatStim(contra);
% contraStim1 = find(contraStimAll == 2);
% contraStim = trialMat1(contra(contraStim1));
stimTrials1 = trialMatStim == 2;
stimTrials = trialMat1(stimTrials1);

noStimTrials1 = trialMatStim == 0;
noStimTrials = trialMat1(noStimTrials1)

%% Assess performance for the highest lums (=255): Need fraction contra choice
highLmvLum = find(dirLums == lumVal_L1);                                   
numelHighLums_L1 = numel(highLmvLum);
corrLhigh1 = trialMat4(highLmvLum);    %correct L-mv trials  
corrLhigh = corrLhigh1 > 1;
corrLhighInds = find(corrLhigh == 1);
numelCorrLhighInds = numel(corrLhighInds);                                  %number correct L trials for 255 lums
fracCorrHighLums_L1 = numelCorrLhighInds/numelHighLums_L1;                  %total fraction correct for 255 lums
numelIncorrLhighInds = numelHighLums_L1 - numelCorrLhighInds;               %number of went-R trials for 255 lums
fracIncorrHighLums_L1 = numelIncorrLhighInds/numelHighLums_L1;              %fraction went-R trials(incorrects)

%255 stimulation trials
highLmvLumStim = trialMatStim(highLmvLum);                                  
numelStimHighL = numel(find(highLmvLumStim == 2));
highLmvStim = highLmvLum(highLmvLumStim==2);                                %these are the stim, 255 lum trials that were correct

corrL1High1 = trialMat4(highLmvStim);                                       %find the ones that were correct
% corr_L1high1 = trialMat4(highLmvLum);                                       %L-rewarded T-hold crossings = correct
%   for i = 1:length(highLmvStim)
%       if highLmvStim(i) > length(corr_L1high1)
%           break
%       else
%          corrhighLmvStim(i) = corr_L1high1(highLmvStim(i));                 %Thold crossings for correct, stimulated, L, high lum trials
%       end
%   end
corrHighLmvStim = isfinite(corrL1High1);
numel_corrHighLmvStim = numel(corrHighLmvStim );
numelcorrHighLmvStim = numel(find(corrHighLmvStim == 1));                   % here is num_L for 255
numelincorrHighLmvStim = numel_corrHighLmvStim - numelcorrHighLmvStim;      % here is num_R (incorrects) for 255 luminance
fracCorr_L1Stim = numelcorrHighLmvStim/numel_corrHighLmvStim;               %This is what I want for correct, stim, high lum % corrects
fracIncorrHighLums_L1Stim = numelincorrHighLmvStim/numelHighLmvStim;   %fraction went-R trials(incorrects)

% numelHighLmvStim = numel(corrhighLmvStim);
% corrHighLmvStim1 = isfinite(corrhighLmvStim);
% numelcorrHighLmvStim = numel(find(corrHighLmvStim1 == 1));             % here is num_L for 255 luminance
% numelincorrHighLmvStim = numelHighLmvStim - numelcorrHighLmvStim;      % here is num_R (incorrects) for 255 luminance
% fracCorr_L1Stim = numelcorrHighLmvStim/numelHighLmvStim;               %This is what I want for correct, stim, high lum % corrects
% fracIncorrHighLums_L1Stim = numelincorrHighLmvStim/numelHighLmvStim;   %fraction went-R trials(incorrects)

%Unstim high lum L trials:
numelUnstimHighL = numel(find(highLmvLumStim == 0));
highLmvUnstim = highLmvLum(highLmvLumStim == 0);                        %these are the stim high lum L trials
%     corr_L1high1 = trialMat4(highLmvLum);
  for i = 1:length(highLmvUnstim)
      if highLmvUnstim(i) > length(corr_L1high1)
          break
      else
         corrhighLmvUnstim(i) = corr_L1high1(highLmvUnstim(i));         %Thold crossings for correct, stimulated, high lum trials
      end
  end
numelHighLmvUnstim = numel(corrhighLmvUnstim);
corrHighLmvUnstim1 = isfinite(corrhighLmvUnstim);
numelcorrHighLmvUnstim = numel(find(corrHighLmvUnstim1 == 1));          % here is num_L
numelincorrHighLmvUnstim = numelHighLmvUnstim - numelcorrHighLmvUnstim; % here is num_R (incorrects) for 255 luminance
fracCorr_L1Unstim = numelcorrHighLmvUnstim/numelHighLmvUnstim;          %This is what I want for correct, unstim, high lum % corrects
fracIncorrHighLums_L1Unstim = numelincorrHighLmvUnstim/numelHighLmvUnstim;

%% all should-R-mv high lum trials (= -255)
highRmvLum = find(dirLums == lumVal_R1);
numelHighLums_R1 = numel(highRmvLum);
corrRhigh1 = trialMat6(highRmvLum);     %correct R-mv trials
corrRhigh = corrRhigh1 > 1;
corrRhighInds = find(corrRhigh == 1);
numelCorr_RhighInds = numel(corrRhighInds);                                 %number correct R trials for -255 lums
fracCorr_R1 = numelCorr_RhighInds/numelHighLums_R1;                         %fraction correct R trials for -255 lums
numelIncorrRhighInds = numelHighLums_R1 - numelCorr_RhighInds;              %number of went-L trials for -255 lums
fracIncorr_R1 = numelIncorrRhighInds/numelHighLums_R1;                      %fraction went-L trials (incorrects)

%separate out stimulated R high lum trials:
highRmvLumStim = trialMatStim(highRmvLum);
numelStimHighR = numel(find(highRmvLumStim == 2));
highRmvStim = highRmvLum(highRmvLumStim==2);                               %these are the stim high lum L trials
for i = 1:length(highRmvStim)
  if highRmvStim(i) > length(corrRhigh1)
      break
  else
     corrhighRmvStim(i) = corrRhigh1(highRmvStim(i));                       %Thold crossings for correct, stimulated, high lum trials
  end
end
numelHighRmvStim = numel(corrhighRmvStim);
corrHighRmvStim1 = isfinite(corrhighRmvStim);
numelcorrHighRmvStim = numel(find(corrHighRmvStim1 == 1));                  %this is num_R for -255
numelincorrHighRmvStim = numelHighRmvStim - numelcorrHighRmvStim;           %this is num_L (incorrects) for -255
fracCorr_R1Stim = numelcorrHighRmvStim/numelHighRmvStim;                    %This is what I want for correct, stim, high lum % corrects
fracIncorr_R1stim = numelincorrHighRmvStim/numelHighRmvStim;                %fraction went_L trials (incorrects)

%Unstim high lum R trials (-255):
numelUnstimHighR = numel(find(highRmvLumStim == 0));
highRmvUnstim = highRmvLum(highRmvLumStim == 0);                            %these are the stim high lum R trials
corr_R1high1 = trialMat6(highRmvLum);
corrhighRmvUnstim = [];
  for i = 1:length(highRmvUnstim)
      if highRmvUnstim(i) > length(corr_R1high1)
          break
      else
         corrhighRmvUnstim(i) = corr_R1high1(highRmvUnstim(i));             %Thold crossings for correct, stimulated, high lum trials
      end
  end
numelHighRmvUnstim = numel(corrhighRmvUnstim);
corrHighRmvUnstim1 = isfinite(corrhighRmvUnstim);
numelcorrHighRmvUnstim = numel(find(corrHighRmvUnstim1 == 1));              %this is num_R for -255 unstim
numelincorrHighRmvUnstim = numelHighRmvUnstim - numelcorrHighRmvUnstim;     %this is num_L for -255 unstim
fracCorr_R1Unstim = numelcorrHighRmvUnstim/numelHighRmvUnstim;              %This is what I want for correct, unstim, high lum R trials' % corrects

%% Assess performance for the mid low lum (= 64):
lowLmvLum = find(dirLums == lumVal_L2);   %trials of lum vals = 65
numelLowLums = numel(lowLmvLum);
corr_L2low1 = trialMat4(lowLmvLum);
corrLlow = corr_L2low1 > 1;
corrLlowInds = find(corrLlow == 1);
numelCorrLlowInds = numel(corrLlowInds);                                    %number correct L trials for 64 lums
fracCorrLowLLums = numelCorrLlowInds/numelLowLums;                          %fraction correct L trials for 64
numelIncorrLlowInds = numelLowLums - numelCorrLlowInds;                     %number of went-R (incorrect) trials for 64 lums
fracIncorrLowLlums = numelIncorrLlowInds/numelLowLums;                      %fraction of went-R (incorrect) trials

% Separate out stim L mid-lum trials
lowLmvLumStim = trialMatStim(lowLmvLum);
numelStimLowL = numel(find(lowLmvLumStim == 2));
lowLmvStim = lowLmvLum(lowLmvLumStim==2);                                   %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
corrlowLmvStim = [];
  for i = 1:length(lowLmvStim)
      if lowLmvStim(i) > length(corr_L2low1)
          break
      else
         corrlowLmvStim(i) = corr_L2low1(lowLmvStim(i));                    %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowLmvStim = numel(corrlowLmvStim);
corrLowLmvStim1 = isfinite(corrlowLmvStim);                                 %total trials of should LowLmvStim
numelcorrLowLmvStim = numel(find(corrLowLmvStim1 == 1));                    %here is num_L for lumValL2
numelincorrLowLmvStim = numelLowLmvStim - numelcorrLowLmvStim;              %here is num_R (incorrects) for lumValL2    
fracCorr_L2Stim = numelcorrLowLmvStim/numelLowLmvStim;                      %This is what I want for correct, stim, low lum % corrects
fracIncorr_L2Stim = numelincorrLowLmvStim/numelLowLmvStim                   %fraction of went-R (incorrect) trials for lum = 64

% Separate out Unstim +64 lum trials:
%     lowLmvLumStim = trialMatStim(lowLmvLum);
numelStimLowL = numel(find(lowLmvLumStim == 0));
lowLmvUnstim = lowLmvLum(lowLmvLumStim == 0);                               %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
corrlowLmvUnstim = [];
  for i = 1:length(lowLmvUnstim)
      if lowLmvUnstim(i) > length(corr_L2low1)
          break
      else
         corrlowLmvUnstim(i) = corr_L2low1(lowLmvUnstim(i));                %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowLmvUnstim = numel(corrlowLmvUnstim);
corrLowLmvUnstim1 = isfinite(corrlowLmvUnstim);
numelcorrLowLmvUnstim = numel(find(corrLowLmvUnstim1 == 1));
numelincorrLowLmvUnstim = numelLowLmvUnstim - numelcorrLowLmvUnstim;        %this is num_R (incorrects) for 64 Lum unstim
fracCorr_L2Unstim = numelcorrLowLmvUnstim/numelLowLmvUnstim;                %this is what I want for correct, stim, low lum % corrects
fracIncorr_L2Unstim = numelincorrLowLmvUnstim/numelLowLmvUnstim;            %fraction of went-R (incorrect) trials for lum = 64

%% all should-R-mv low lum trials (= -64)
lowRmvLum = find(dirLums == lumVal_R2);
numelLowR_inds = numel(lowRmvLum);
corr_R2low1 = trialMat6(lowRmvLum);
corr_Rlow = corr_R2low1 > 1;
corr_R2lowInds = find(corr_Rlow == 1);
numelCorr_R2lowInds = numel(corr_R2lowInds);                                %number correct R trials for -64 lums
fracCorr_R2 = numelCorr_R2lowInds/numelLowR_inds;                           %fraction correct L low lums
numelIncorrRlowInds = numelLowR_inds - numelCorr_R2lowInds;                 %number of went-L (incorrect) trials for -64 lums
fracIncorr_R2 = numelIncorrRlowInds/numelLowR_inds;                         %frac of went-L (incorrect) trials

%separate out stimulated R mid lum (-64) trials:
lowRmvLumStim = trialMatStim(lowRmvLum);
numelStimLowR = numel(find(lowRmvLumStim == 2));
lowRmvStim = lowRmvLum(lowRmvLumStim==2);                                   %these are the stim  low lum L trials
corr_R2low1 = trialMat6(lowRmvLum);
corrlowRmvStim = [];
for i = 1:length(lowRmvStim)
  if lowRmvStim(i) > length(corr_R2low1)
      break
  else
     corrlowRmvStim(i) = corr_R2low1(lowRmvStim(i));                        %Thold crossings for correct, stimulated, low lum trials
  end
end
numelLowRmvStim = numel(corrlowRmvStim);
corrLowRmvStim1 = isfinite(corrlowRmvStim);
numelcorrLowRmvStim = numel(find(corrLowRmvStim1 == 1));                    %this is the num_R 
numelincorrLowRmvStim = numelLowRmvStim - numelcorrLowRmvStim;              %this is the num_L (incorrects) for 64 lum on stim trials
fracCorr_R2Stim = numelcorrLowRmvStim/numelLowRmvStim;                      %This is what I want for correct, stim, low lum % corrects

%Unstim mid lum R (-64) trials:
numelUnstimLowR = numel(find(lowRmvLumStim == 0));
lowRmvUnstim = lowRmvLum(lowRmvLumStim == 0);                               %these are the stim low lum R trials
corr_R2low1 = trialMat6(lowRmvLum);
corrlowRmvUnstim = [];
  for i = 1:length(lowRmvUnstim)
      if lowRmvUnstim(i) > length(corr_R2low1)
          break
      else
         corrlowRmvUnstim(i) = corr_R2low1(lowRmvUnstim(i));                %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowRmvUnstim = numel(corrlowRmvUnstim);
corrLowRmvUnstim1 = isfinite(corrlowRmvUnstim);
numelcorrLowRmvUnstim = numel(find(corrLowRmvUnstim1 == 1));                %this is the num_R
numelincorrLowRmvUnstim = numelLowRmvUnstim - numelcorrLowRmvUnstim;        %this is num_L (incorrects) for unstim -64 lum
fracCorr_R2Unstim = numelcorrLowRmvUnstim/numelLowRmvUnstim;                %This is what I want for correct, unstim, low lum R trials' % corrects

%% Assess performance for the 50/50 lowest lum (=15):
lowestLmvLum = find(dirLums == lumVal_L3);  %trials of lum vals = 15
numelLowestLLums = numel(lowestLmvLum);
corrLlowest1 = trialMat4(lowestLmvLum );
corrLlowest = corrLlowest1 > 1;
corrLlowestInds = find(corrLlowest == 1);
numelCorrLlowestInds = numel(corrLlowestInds);                              %number correct L trials for 15 lums
fracCorrLowestLLums = numelCorrLlowestInds/numelLowestLLums;                 %fraction correct L trials for 15 lums
numelIncorrLlowestInds = numelLowestLLums - numelCorrLlowestInds;           %number of incorrect (went-R) trials for 15 lums
fracIncorrLowestLlums = numelIncorrLlowestInds/numelLowestLLums;            %frac of incorrect (went-R) trials

% Separate out stim 15 lum trials
lowestLmvLumStim = trialMatStim(lowestLmvLum);
numelStimLowestL = numel(find(lowestLmvLumStim == 2));
lowestLmvStim = lowestLmvLum(lowestLmvLumStim==2);                          %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
corrlowestLmvStim = [];
  for i = 1:length(lowestLmvStim)
      if lowestLmvStim(i) > length(corrLlowest1)
          break
      else
         corrlowestLmvStim(i) = corrLlowest1(lowestLmvStim(i));            %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowestLmvStim = numel(corrlowestLmvStim);
corrLowestLmvStim1 = isfinite(corrlowestLmvStim);                           %total trials of should LowLmvStim
numelcorrLowestLmvStim = numel(find(corrLowestLmvStim1 == 1));              %here is num_L for lumValL2
numelincorrLowestLmvStim = numelLowestLmvStim - numelcorrLowestLmvStim;     %here is num_R (incorrects) for lumValL2    
fracCorr_L3Stim = numelcorrLowestLmvStim/numelLowestLmvStim;                %This is what I want for correct, stim, low lum % corrects
fracIncorr_L3Stim = numelincorrLowestLmvStim/numelLowestLmvStim;            %fraction of went-R (incorrect) trials for lum = 64

% Separate out unstim 15 lum trials
lowestLmvLumUnstim = trialMatStim(lowestLmvLum);
numelUnstimLowestL = numel(find(lowestLmvLumUnstim == 2));
lowestLmvUnstim = lowestLmvLum(lowestLmvLumUnstim==2);                     %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
corrlowestLmvUnstim = [];
  for i = 1:length(lowestLmvUnstim)
      if lowestLmvUnstim(i) > length(corrLlowest1)
          break
      else
         corrlowestLmvUnstim(i) = corrLlowest1(lowestLmvUnstim(i));         %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowestLmvUnstim = numel(corrlowestLmvUnstim);
corrLowestLmvUnstim1 = isfinite(corrlowestLmvUnstim);                       %total trials of should LowLmvStim
numelcorrLowestLmvUnstim = numel(find(corrLowestLmvUnstim1 == 1));          %here is num_L for lumValL2
numelincorrLowestLmvUnstim = numelLowestLmvUnstim - numelcorrLowestLmvUnstim; %here is num_R (incorrects) for lumValL2    
fracCorr_L3Unstim = numelcorrLowestLmvUnstim/numelLowestLmvUnstim;          %This is what I want for correct, stim, low lum % corrects
fracIncorr_L3Unstim = numelincorrLowestLmvUnstim/numelLowestLmvUnstim;      %fraction of went-R (incorrect) trials for lum = 64

%% Assess performance for the 50/50 lowest lum (= -15):
lowestRmvLum = find(dirLums == lumVal_R3);
numelLowestR_inds = numel(lowestRmvLum);
corr_Rlowest1 = trialMat6(lowestRmvLum);
corr_Rlowest = corr_Rlowest1 > 1;
corr_RlowestInds = find(corr_Rlowest == 1);
numelCorr_RlowestInds = numel(corr_RlowestInds);                            %number correct R trials for -15 lums
fracCorrLowestRLums = numelCorr_RlowestInds/numelLowestR_inds;              %fraction correct R lowest lums (-15)
numelIncorrRlowestInds = numelLowestR_inds - numelCorr_RlowestInds;         %number went-L (incorrect) trials
fracIncorrRlowestLums = numelIncorrRlowestInds/numelLowestR_inds;           %frac went-L (incorrect) trials

%separate out stimulated R (-15) lowest lum trials:
lowestRmvLumStim = trialMatStim(lowestRmvLum);
numelStimLowestR = numel(find(lowestRmvLumStim == 2));
lowestRmvStim = lowestRmvLum(lowestRmvLumStim==2);                          %these are the stim  lowrst lum L trials
corr_R3lowest1 = trialMat6(lowestRmvLum);
corrlowestRmvStim = [];
for i = 1:length(lowestRmvStim)
  if lowestRmvStim(i) > length(corr_R3lowest1)
      break
  else
     corrlowestRmvStim(i) = corr_R3lowest1(lowestRmvStim(i));               %Thold crossings for correct, stimulated, lowest lum trials
  end
end
numelLowestRmvStim = numel(corrlowestRmvStim);
corrLowestRmvStim1 = isfinite(corrlowestRmvStim);
numelcorrLowestRmvStim = numel(find(corrLowestRmvStim1 == 1));              %num_R for -15 stim
numelincorrLowestRmvStim =  numelLowestRmvStim - numelcorrLowestRmvStim;    %num_L for -15 stim
fracCorr_R3Stim = numelcorrLowestRmvStim/numelLowestRmvStim;                %This is what I want for correct, stim, lowest lum % corrects

%Unstim lowest lum R (-15) trials:
numelUnstimLowR = numel(find(lowestRmvLumStim == 0));
lowestRmvUnstim = lowestRmvLum(lowestRmvLumStim == 0);                      %these are the stim lowest lum R trials
corr_R3lowest1 = trialMat6(lowestRmvLum);
corrlowestRmvUnstim = [];
  for i = 1:length(lowestRmvUnstim)
      if lowestRmvUnstim(i) > length(corr_R3lowest1)
          break
      else
         corrlowestRmvUnstim(i) = corr_R3lowest1(lowestRmvUnstim(i));       %Thold crossings for correct, stimulated, low lum trials
      end
  end
numelLowestRmvUnstim = numel(corrlowestRmvUnstim);
corrLowestRmvUnstim1 = isfinite(corrlowestRmvUnstim);
numelcorrLowestRmvUnstim = numel(find(corrLowestRmvUnstim1 == 1));          %num_R for -15 unstim
numelincorrLowestRmvUnstim = numelLowestRmvUnstim - numelcorrLowestRmvUnstim;%num_L for -15 unstim
fracCorr_R3Unstim = numelcorrLowestRmvUnstim/numelLowestRmvUnstim;          %This is what I want for correct, unstim, low lum R trials' % corrects


  %% Plot % correct 
  %for the total plot:  works for 8.23.17 data
    figure; hold on;
    xAll = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L3, lumVal_L2, lumVal_L1]'; 
    yAll = [fracCorr_R1, fracCorr_R2, fracCorrLowestRLums, fracIncorrLowestLlums, fracIncorrLowLlums, fracIncorrHighLums_L1]';
    plotAll = plot(xAll, yAll, 'k.', 'MarkerSize', 18);
        set(gca, 'XLim', [-255 255]);
%         set(gca, 'YLim', [0 1]);
        xlabel('Luminance vals', 'FontSize', 14);
        ylabel('Fraction contra choice', 'FontSize', 14);
     num_L_All = [numelIncorrRhighInds, numelIncorrRlowInds, numelIncorrRlowestInds, numelCorrLlowestInds, numelCorrLlowInds, numelCorrLhighInds]';
     num_R_All = [numelCorr_RhighInds, numelCorr_R2lowInds, numelCorr_RlowestInds, numelIncorrLlowestInds, numelIncorrLlowInds, numelIncorrLhighInds]';
    
    best_fit_All = glmfit(xAll, [num_R_All (num_R_All + num_L_All)], 'binomial');
    x_axis_fit = -255:255;
    y_fit = glmval(best_fit_All, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fit, 'k-');
    
    mouseInfo_info = strcat(mouseName, ',', num2str(numSession), 'allTrials');
    title(mouseInfo_info);
    
    %% for the stim v unstim plot:    
    figure; hold on;
    x2 = [-257 -66 -17 17 66 257]; %to jitter the unstim results so they don't plot over stim
    yUnstim = [fracCorr_R1Unstim, fracCorr_R2Unstim, fracCorr_R3Unstim, fracIncorr_L3Unstim, fracIncorr_L2Unstim, fracIncorrHighLums_L1Unstim]';   %For y axis = % correct R 
    plotUnstim = plot(x2, yUnstim, 'k.', 'MarkerSize', 18);
    hold on;
    num_L_unstim = [numelincorrHighRmvUnstim, numelincorrLowRmvUnstim, numelincorrLowestRmvUnstim, numelcorrLowestLmvUnstim, numelcorrLowLmvUnstim, numelcorrHighLmvUnstim]';
    num_R_unstim = [numelcorrHighRmvUnstim, numelcorrLowRmvUnstim, numelcorrLowestRmvUnstim, numelincorrLowestLmvUnstim, numelincorrLowLmvUnstim, numelincorrHighLmvUnstim]';
    best_fit_Unstim = glmfit(x2, [num_R_unstim (num_R_unstim + num_L_unstim)], 'binomial');
    x_axis_fit = -255:255;
    y_fit = glmval(best_fit_Unstim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fit, 'k-');

    %Stim trials
    xStim = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L3, lumVal_L2, lumVal_L1]'; 
    yStim = [fracCorr_R1Stim, fracCorr_R2Stim, fracCorr_R3Stim, fracIncorr_L3Stim, fracIncorr_L2Stim, fracIncorrHighLums_L1Stim]';   
    plotStim = plot(xStim, yStim, 'bl.', 'MarkerSize', 18);

    set(gca, 'XLim', [-260 260], 'Ylim', [0 1]);
        xlabel('Luminance vals', 'FontSize', 14);
        ylabel('Fraction correct contra mv', 'FontSize', 14);
    num_L_stim = [numelincorrHighLmvStim, numelincorrLowRmvStim, numelincorrLowestRmvStim, numelcorrLowestLmvStim, numelcorrLowLmvStim, numelcorrHighLmvStim]';
    num_R_stim = [numelcorrHighRmvStim, numelcorrLowRmvStim, numelcorrLowestRmvStim, numelincorrLowestLmvStim, numelincorrLowLmvStim, numelincorrHighLmvStim]';

    best_fit_stim = glmfit(xStim, [num_R_stim (num_R_stim + num_L_stim)], 'binomial');
    x_axis_fit = -255:255;
    y_fitStim = glmval(best_fit_stim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fitStim, 'bl-');

    mouseInfo = strcat(mouseName, ',', num2str(numSession));
    title(mouseInfo);

