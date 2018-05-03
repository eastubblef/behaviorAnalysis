
function [alignTimes, positions] = lums(taskbase, contra, ipsi, lumTrials)

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
numSession = 170823;

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

%% Assess performance for the highest lums (=255):
highLums = find(lumTrials == lumVal_L1);                                    %note these are trials that are either -255 or 255
numelHighLums = numel(highLums);
corrLhigh1 = trialMat4(highLums);  %correct L trials  
corrRhigh1 = trialMat6(highLums);  %correct R trials
corrLhigh = corrLhigh1 > 1;
corrLhighInds = find(corrLhigh == 1);
numelCorrLhighInds = numel(corrLhighInds);
corrRhigh = corrRhigh1 > 1;
corrRhighInds = find(corrRhigh == 1);
numelCorrRhighInds = numel(corrRhighInds);
numelCorrBothHigh = numelCorrLhighInds + numelCorrRhighInds;
fracCorrHighLums = numelCorrBothHigh/numelHighLums;                         %total fraction correct for both contra/ipsi trials

%separate out stim L high lum trials:
highLmvLum = find(dirLums == lumVal_L1);
numelHighL_inds = numel(highLmvLum);
    highLmvLumStim = trialMatStim(highLmvLum);
    numelStimHighL = numel(find(highLmvLumStim == 2));
    highLmvStim = highLmvLum(highLmvLumStim==2);                            %these are the stim high lum L trials
corr_L1high1 = trialMat4(highLmvLum);
      for i = 1:length(highLmvStim)
          if highLmvStim(i) > length(corr_L1high1)
              break
          else
             corrhighLmvStim(i) = corr_L1high1(highLmvStim(i));             %Thold crossings for correct, stimulated, L, high lum trials
          end
      end
    numelHighLmvStim = numel(corrhighLmvStim);
    corrHighLmvStim1 = isfinite(corrhighLmvStim);
    numelcorrHighLmvStim = numel(find(corrHighLmvStim1 == 1));             % here is num_L for 255 luminance
    numelincorrHighLmvStim = numelHighLmvStim - numelcorrHighLmvStim;      % here is num_R (incorrects) for 255 luminance
    fracCorr_L1Stim = numelcorrHighLmvStim/numelHighLmvStim;               %This is what I want for correct, stim, high lum % corrects
    plot_L1Stim = 1 - fracCorr_L1Stim;   %For y axis = % correct R 

%All L high lum trials
corr_Lhigh = corr_L1high1 > 1;
corr_L1highInds = find(corr_Lhigh == 1);
numelCorr_L1highInds = numel(corr_L1highInds);
fracCorr_L1 = numelCorr_L1highInds/numelHighL_inds;                         %fraction correct L highest lums

%Unstim high lum L trials:
%     highLmvLumStim = trialMatStim(highLmvLum);
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

    fracCorr_L1Unstim = numelcorrHighLmvUnstim/numelHighLmvUnstim;         %This is what I want for correct, unstim, high lum % corrects
    plot_L1Unstim = 1 - fracCorr_L1Unstim;   %For y axis = % correct R 

%all R mv high lum trials
highRmvLum = find(dirLums == lumVal_R1);
numelHighR_inds = numel(highRmvLum);
corr_R1high1 = trialMat6(highRmvLum);
corr_Rhigh = corr_R1high1 > 1;
corr_R1highInds = find(corr_Rhigh == 1);
numelCorr_R1highInds = numel(corr_R1highInds);
fracCorr_R1 = numelCorr_R1highInds/numelHighR_inds;                         %fraction correct L highest lums

%separate out stimulated R high lum trials:
    highRmvLumStim = trialMatStim(highRmvLum);
    numelStimHighR = numel(find(highRmvLumStim == 2));
    highRmvStim = highRmvLum(highRmvLumStim==2);                            %these are the stim high lum L trials
    for i = 1:length(highRmvStim)
      if highRmvStim(i) > length(corr_R1high1)
          break
      else
         corrhighRmvStim(i) = corr_R1high1(highRmvStim(i));                 %Thold crossings for correct, stimulated, high lum trials
      end
    end
    numelHighRmvStim = numel(corrhighRmvStim);
    corrHighRmvStim1 = isfinite(corrhighRmvStim);
    numelcorrHighRmvStim = numel(find(corrHighRmvStim1 == 1));              %this is num_R for -255
    numelincorrHighRmvStim = numelHighRmvStim - numelcorrHighRmvStim;       %this is num_L (incorrects) for -255
    fracCorr_R1Stim = numelcorrHighRmvStim/numelHighRmvStim;                %This is what I want for correct, stim, high lum % corrects

    %Unstim high lum R trials:
    numelUnstimHighR = numel(find(highRmvLumStim == 0));
    highRmvUnstim = highRmvLum(highRmvLumStim == 0);                        %these are the stim high lum R trials
    corr_R1high1 = trialMat6(highRmvLum);
    corrhighRmvUnstim = [];
      for i = 1:length(highRmvUnstim)
          if highRmvUnstim(i) > length(corr_R1high1)
              break
          else
             corrhighRmvUnstim(i) = corr_R1high1(highRmvUnstim(i));         %Thold crossings for correct, stimulated, high lum trials
          end
      end
    numelHighRmvUnstim = numel(corrhighRmvUnstim);
    corrHighRmvUnstim1 = isfinite(corrhighRmvUnstim);
    numelcorrHighRmvUnstim = numel(find(corrHighRmvUnstim1 == 1));         %this is num_R for -255 unstim
    numelincorrHighRmvUnstim = numelHighRmvUnstim - numelcorrHighRmvUnstim;%this is num_L for -255 unstim
    fracCorr_R1Unstim = numelcorrHighRmvUnstim/numelHighRmvUnstim;         %This is what I want for correct, unstim, high lum R trials' % corrects

% Assess performance for the mid low lum (=65):
lowLums = find(lumTrials == lumVal_L2);   %trials of lum vals = 65
numelLowLums = numel(lowLums);
corrLlow1 = trialMat4(lowLums);
corrRlow1 = trialMat6(lowLums);
corrLlow = corrLlow1 > 1;
corrLlowInds = find(corrLlow == 1);
numelCorrLlowInds = numel(corrLlowInds);
corrRlow = corrRlow1 > 1;
corrRlowInds = find(corrRlow == 1);
numelCorrRlowInds = numel(corrRlowInds);
numelCorrBothlow = numelCorrLlowInds + numelCorrRlowInds
fracCorrLowLums = numelCorrBothlow/numelLowLums;                            %total fraction correct for both contra/ipsi trials

lowLmvLum = find(dirLums == lumVal_L2);
numelLowL_inds = numel(lowLmvLum);
corr_L2low1 = trialMat4(lowLmvLum);
corr_Llow = corr_L2low1 > 1;
corr_L2lowInds = find(corr_Llow == 1);
numelCorr_L2lowInds = numel(corr_L2lowInds);
fracCorr_L2 = numelCorr_L2lowInds/numelLowL_inds;                           %fraction correct L low lums

lowRmvLum = find(dirLums == lumVal_R2);
numelLowR_inds = numel(lowRmvLum);
corr_R2low1 = trialMat6(lowRmvLum);
corr_Rlow = corr_R2low1 > 1;
corr_R2lowInds = find(corr_Rlow == 1);
numelCorr_R2lowInds = numel(corr_R2lowInds);
fracCorr_R2 = numelCorr_R2lowInds/numelLowR_inds;                           %fraction correct L low lums

% Separate out stim L mid-lum trials
    lowLmvLumStim = trialMatStim(lowLmvLum);
    numelStimLowL = numel(find(lowLmvLumStim == 2));
    lowLmvStim = lowLmvLum(lowLmvLumStim==2);                               %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
    corrlowLmvStim = [];
      for i = 1:length(lowLmvStim)
          if lowLmvStim(i) > length(corr_L2low1)
              break
          else
             corrlowLmvStim(i) = corr_L2low1(lowLmvStim(i));                %Thold crossings for correct, stimulated, low lum trials
          end
      end
    numelLowLmvStim = numel(corrlowLmvStim);
    corrLowLmvStim1 = isfinite(corrlowLmvStim);                             %total trials of should LowLmvStim
    numelcorrLowLmvStim = numel(find(corrLowLmvStim1 == 1));                %here is num_L for lumValL2
    numelincorrLowLmvStim = numelLowLmvStim - numelcorrLowLmvStim;          %here is num_R (incorrects) for lumValL2
    
    fracCorr_L2Stim = numelcorrLowLmvStim/numelLowLmvStim;                  %This is what I want for correct, stim, low lum % corrects
    plot_L2Stim = 1 - fracCorr_L2Stim;                                      %For y axis = % correct R 

% Separate out Unstim L mid-lum trials
%     lowLmvLumStim = trialMatStim(lowLmvLum);
    numelStimLowL = numel(find(lowLmvLumStim == 0));
    lowLmvUnstim = lowLmvLum(lowLmvLumStim == 0);                           %these are the stim low lum L trials
%     corr_L2low1 = trialMat4(lowLmvLum);
    corrlowLmvUnstim = [];
      for i = 1:length(lowLmvUnstim)
          if lowLmvUnstim(i) > length(corr_L2low1)
              break
          else
             corrlowLmvUnstim(i) = corr_L2low1(lowLmvUnstim(i));            %Thold crossings for correct, stimulated, low lum trials
          end
      end
          
    numelLowLmvUnstim = numel(corrlowLmvUnstim);
    corrLowLmvUnstim1 = isfinite(corrlowLmvUnstim);
    numelcorrLowLmvUnstim = numel(find(corrLowLmvUnstim1 == 1));
    numelincorrLowLmvUnstim = numelLowLmvUnstim - numelcorrLowLmvUnstim;   %this is num_R (incorrects) for 64 Lum unstim
    fracCorr_L2Unstim = numelcorrLowLmvUnstim/numelLowLmvUnstim;           %this is what I want for correct, stim, low lum % corrects
    plot_L2Unstim = 1 - fracCorr_L2Unstim;   %For y axis = % correct R 

%separate out stimulated R mid lum trials:
    lowRmvLumStim = trialMatStim(lowRmvLum);
    numelStimLowR = numel(find(lowRmvLumStim == 2));
    lowRmvStim = lowRmvLum(lowRmvLumStim==2);                               %these are the stim  low lum L trials
    corr_R2low1 = trialMat6(lowRmvLum);
    corrlowRmvStim = [];
    for i = 1:length(lowRmvStim)
      if lowRmvStim(i) > length(corr_R2low1)
          break
      else
         corrlowRmvStim(i) = corr_R2low1(lowRmvStim(i));                    %Thold crossings for correct, stimulated, low lum trials
      end
    end
    numelLowRmvStim = numel(corrlowRmvStim);
    corrLowRmvStim1 = isfinite(corrlowRmvStim);
    numelcorrLowRmvStim = numel(find(corrLowRmvStim1 == 1));                %this is the num_L 
    numelincorrLowRmvStim = numelLowRmvStim - numelcorrLowRmvStim;          %this is the num_R (incorrects) for 64 lum on stim trials
    fracCorr_R2Stim = numelcorrLowRmvStim/numelLowRmvStim;                  %This is what I want for correct, stim, low lum % corrects

    %Unstim mid lum R trials:
    numelUnstimLowR = numel(find(lowRmvLumStim == 0));
    lowRmvUnstim = lowRmvLum(lowRmvLumStim == 0);                           %these are the stim low lum R trials
    corr_R2low1 = trialMat6(lowRmvLum);
    corrlowRmvUnstim = [];
      for i = 1:length(lowRmvUnstim)
          if lowRmvUnstim(i) > length(corr_R2low1)
              break
          else
             corrlowRmvUnstim(i) = corr_R2low1(lowRmvUnstim(i));            %Thold crossings for correct, stimulated, low lum trials
          end
      end
    numelLowRmvUnstim = numel(corrlowRmvUnstim);
    corrLowRmvUnstim1 = isfinite(corrlowRmvUnstim);
    numelcorrLowRmvUnstim = numel(find(corrLowRmvUnstim1 == 1));
    numelincorrLowRmvUnstim = numelLowRmvUnstim - numelcorrLowRmvUnstim;   %this is num_L (incorrects) for unstim 64 lum
    fracCorr_R2Unstim = numelcorrLowRmvUnstim/numelLowRmvUnstim;           %This is what I want for correct, unstim, low lum R trials' % corrects

%% Assess performance for the 50/50 lowest lum (=15):
lowestLums = find(lumTrials == lumVal_L3);  %trials of lum vals = 15
numelLowestLums = numel(lowestLums);
corrLlowest1 = trialMat4(lowestLums);
corrRlowest1 = trialMat6(lowestLums);
corrLlowest = corrLlowest1 > 1;
corrLlowestInds = find(corrLlowest == 1);
numelCorrLlowestInds = numel(corrLlowestInds);
corrRlowest = corrRlowest1 > 1;
corrRlowestInds = find(corrRlowest == 1);
numelCorrRlowestInds = numel(corrRlowestInds);
numelCorrBothlowest = numelCorrLlowestInds + numelCorrRlowestInds;
fracCorrLowestLums = numelCorrBothlowest/numelLowestLums;                  %total fraction correct for both contra/ipsi trials

lowestLmvLum = find(dirLums == lumVal_L3);
numelLowestL_inds = numel(lowestLmvLum);
corr_L3lowest1 = trialMat4(lowestLmvLum);
corr_Llowest = corr_L3lowest1 > 1;
corr_L3lowestInds = find(corr_Llowest == 1);
numelCorr_L3lowestInds = numel(corr_L3lowestInds);
fracCorr_L3 = numelCorr_L3lowestInds/numelLowestL_inds;                    %fraction correct L hlowest lums

lowestRmvLum = find(dirLums == lumVal_R3);
numelLowestR_inds = numel(lowestRmvLum);
corr_R3lowest1 = trialMat6(lowestRmvLum);
corr_Rlowest = corr_R3lowest1 > 1;
corr_R3lowestInds = find(corr_Rlowest == 1);
numelCorr_R3lowestInds = numel(corr_R3lowestInds);
fracCorr_R3 = numelCorr_R3lowestInds/numelLowestR_inds;                    %fraction correct L lowest lums

% Separate out stim L lowest-lum trials
    lowestLmvLumStim = trialMatStim(lowestLmvLum);
    numelStimLowestL = numel(find(lowestLmvLumStim == 2));
    lowestLmvStim = lowestLmvLum(lowestLmvLumStim==2);                     %these are the stim lowest lum L trials
%    corr_R3lowest1 = trialMat6(lowestRmvLum);
    corrlowestLmvStim = [];
      for i = 1:length(lowestLmvStim)
          if lowestLmvStim(i) > length(corr_R3lowest1)
              break
          else
             corrlowestLmvStim(i) = corr_L3lowest1(lowestLmvStim(i));      %Thold crossings for correct, stimulated, lowest lum trials
          end
      end
    numelLowestLmvStim = numel(corrlowestLmvStim);
    corrLowestLmvStim1 = isfinite(corrlowestLmvStim);
    numelcorrLowestLmvStim = numel(find(corrLowestLmvStim1 == 1));
    numelincorrLowestLmvStim = numelLowestLmvStim - numelcorrLowestLmvStim;
    fracCorr_L3Stim = numelcorrLowestLmvStim/numelLowestLmvStim;           %This is what I want for correct, stim, lowest lum % corrects
    plot_L3Stim = 1 - fracCorr_L3Stim;   %For y axis = % correct R 

% Separate out Unstim L lowest-lum trials
%     lowLmvLumStim = trialMatStim(lowLmvLum);
    numelUnstimLowestL = numel(find(lowestLmvLumStim == 0));
    lowestLmvUnstim = lowestLmvLum(lowestLmvLumStim == 0);                  %these are the stim high lum L trials
%    corr_L3lowest1 = trialMat4(lowestLmvLum);
    corrlowestLmvUnstim = [];
      for i = 1:length(lowLmvUnstim)
          if lowestLmvUnstim(i) > length(corr_L3lowest1)
              break
          else
             corrlowestLmvUnstim(i) = corr_L3lowest1(lowestLmvUnstim(i));  %Thold crossings for correct, stimulated, high lum trials
          end
      end
          
    numelLowestLmvUnstim = numel(corrlowestLmvUnstim);
    corrLowestLmvUnstim1 = isfinite(corrlowestLmvUnstim);
    numelcorrLowestLmvUnstim = numel(find(corrLowestLmvUnstim1 == 1));
    numelincorrLowestLmvUnstim = numelLowestLmvUnstim - numelcorrLowestLmvUnstim;
    fracCorr_L3Unstim = numelcorrLowestLmvUnstim/numelLowestLmvUnstim;     %This is what I want for correct, stim, high lum % corrects
    plot_L3Unstim = 1 - fracCorr_L3Unstim;                                 %For y axis = % correct R 

%separate out stimulated R lowest lum trials:
    lowestRmvLumStim = trialMatStim(lowestRmvLum);
    numelStimLowestR = numel(find(lowestRmvLumStim == 2));
    lowestRmvStim = lowestRmvLum(lowestRmvLumStim==2);                        %these are the stim  lowrst lum L trials
    corr_R3lowest1 = trialMat6(lowestRmvLum);
    corrlowestRmvStim = [];
    for i = 1:length(lowestRmvStim)
      if lowestRmvStim(i) > length(corr_R3lowest1)
          break
      else
         corrlowestRmvStim(i) = corr_R3lowest1(lowestRmvStim(i));          %Thold crossings for correct, stimulated, lowest lum trials
      end
    end
    numelLowestRmvStim = numel(corrlowestRmvStim);
    corrLowestRmvStim1 = isfinite(corrlowestRmvStim);
    numelcorrLowestRmvStim = numel(find(corrLowestRmvStim1 == 1));
    numelincorrLowestRmvStim =  numelLowestRmvStim - numelcorrLowestRmvStim;
    fracCorr_R3Stim = numelcorrLowestRmvStim/numelLowestRmvStim;           %This is what I want for correct, stim, lowest lum % corrects
    
    %Unstim lowest lum R trials:
    numelUnstimLowR = numel(find(lowestRmvLumStim == 0));
    lowestRmvUnstim = lowestRmvLum(lowestRmvLumStim == 0);                 %these are the stim lowest lum R trials
    corr_R3lowest1 = trialMat6(lowestRmvLum);
    corrlowestRmvUnstim = [];
      for i = 1:length(lowestRmvUnstim)
          if lowestRmvUnstim(i) > length(corr_R3lowest1)
              break
          else
             corrlowestRmvUnstim(i) = corr_R3lowest1(lowestRmvUnstim(i));  %Thold crossings for correct, stimulated, low lum trials
          end
      end
    numelLowestRmvUnstim = numel(corrlowestRmvUnstim);
    corrLowestRmvUnstim1 = isfinite(corrlowestRmvUnstim);
    numelcorrLowestRmvUnstim = numel(find(corrLowestRmvUnstim1 == 1));
    numelincorrLowestRmvUnstim = numelLowestRmvUnstim - numelcorrLowestRmvUnstim;
    fracCorr_R3Unstim = numelcorrLowestRmvUnstim/numelLowestRmvUnstim;     %This is what I want for correct, unstim, low lum R trials' % corrects
 
  %% Plot % correct 
  %for the total plot:  works for 8.23.17 data
    figure; hold on;
    xAll = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L2, lumVal_L1]'; 
    yAll = [fracCorr_R1, fracCorr_R2, fracCorr_R3, 1-fracCorr_L2, 1-fracCorr_L1]';
    plotAll = plot(xAll, yAll, 'k.', 'MarkerSize', 18);
        set(gca, 'XLim', [-255 255]);
%         set(gca, 'YLim', [0 1]);
        xlabel('Luminance vals', 'FontSize', 14);
        ylabel('Fraction correct R-mv', 'FontSize', 14);
%      num_L_All = [numelCorrLhighInds, numelCorrLlowInds, numelCorrLlowestInds, numelCorrRlowInds, numelCorrRhighInds]';
%      num_R_All = [numelCorrRhighInds, numelCorrRlowInds, numelCorrRlowestInds, numelCorrLlowInds, numelCorrLhighInds]';
    num_R_All = [9.5, 9, 5.5, 3.5 2]';
    num_L_All = [1.5, 1, 4.5, 6.5 8]';  %works for 8.23.17 data
    best_fit_All = glmfit(xAll, [num_R_All (num_R_All + num_L_All)], 'binomial');
    x_axis_fit = -255:255;
    y_fit = glmval(best_fit_All, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fit, 'k-');
    
    mouseInfo_info = strcat(mouseName, ',', num2str(numSession), 'allTrials');
    title(mouseInfo_info);
    
    %% for the stim v unstim plot:    
    figure; hold on;
    x2 = [-257 -66 -17 17 66 257]; %to jitter the unstim results so they don't plot over stim
    yUnstim = [fracCorr_R1Unstim, fracCorr_R2Unstim, fracCorr_R3Unstim, plot_L3Unstim, plot_L2Unstim, plot_L1Unstim]';   %For y axis = % correct R 
    plotUnstim = plot(x2, yUnstim, 'k.', 'MarkerSize', 18);
    hold on;
    num_L_unstim = [numelincorrHighLmvUnstim, numelincorrLowLmvUnstim,  numelincorrLowestLmvUnstim, numelcorrLowestLmvUnstim, numelcorrLowLmvUnstim, numelcorrHighLmvUnstim]';
    num_R_unstim = [numelcorrHighRmvUnstim, numelcorrLowRmvUnstim, numelcorrLowestRmvUnstim, numelincorrLowestRmvUnstim, numelincorrLowRmvUnstim, numelincorrHighRmvUnstim]';
    best_fit_Unstim = glmfit(x2, [num_R_unstim (num_R_unstim + num_L_unstim)], 'binomial');
    x_axis_fit = -255:255;
    y_fit = glmval(best_fit_Unstim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fit, 'k-');

    %Stim trials
    xStim = [lumVal_R1, lumVal_R2, lumVal_R3, lumVal_L3, lumVal_L2, lumVal_L1]'; 
    yStim = [fracCorr_R1Stim, fracCorr_R2Stim, fracCorr_R3Stim, plot_L3Stim, plot_L2Stim, plot_L1Stim]';   
    plotStim = plot(xStim, yStim, 'bl.', 'MarkerSize', 18);

    set(gca, 'XLim', [-260 260], 'Ylim', [0 1]);
        xlabel('Luminance vals', 'FontSize', 14);
        ylabel('Fraction correct contra mv', 'FontSize', 14);
%         fiftyL = (numelincorrLowestLmvStim + numelcorrLowestLmvStim)/2;
%         fiftyR = (numelcorrLowestRmvStim + numelincorrLowestLmvStim)/2;
    num_L_stim = [numelincorrHighLmvStim, numelincorrLowLmvStim, numelincorrLowestLmvStim, numelcorrLowestLmvStim, numelcorrLowLmvStim, numelcorrHighLmvStim]';
    num_R_stim = [numelcorrHighRmvStim, numelcorrLowRmvStim, numelcorrLowestRmvStim, numelincorrLowestRmvStim, numelincorrLowRmvStim, numelincorrHighRmvStim]';
%     num_L_stim = [numelincorrHighLmvStim, numelincorrLowLmvStim, fiftyL, numelcorrLowLmvStim, numelcorrHighLmvStim]';
%     num_R_stim = [numelcorrHighRmvStim, numelcorrLowRmvStim, fiftyR, numelincorrLowRmvStim, numelincorrHighRmvStim]';

    best_fit_stim = glmfit(xStim, [num_R_stim (num_R_stim + num_L_stim)], 'binomial');
    x_axis_fit = -255:255;
    y_fitStim = glmval(best_fit_stim, x_axis_fit, 'logit');
    plot(x_axis_fit, y_fitStim, 'bl-');

    mouseInfo = strcat(mouseName, ',', num2str(numSession));
    title(mouseInfo);

    %% Deprecated code
% % trialsToPlot1 = 10;
% % trialsToPlot2 = 20; 
% % trialsToPlot3 = 30; 
% % 
% %//////////////////////////////////////////////////////////////////////////
% %fileDirectory = '/Users/stubblefielde/Desktop/VGATone/easiest/moreTrials/'    %for all trials (4/15/15) 
% % cd '/Users/stubblefielde/Desktop/VGATone/easiest/150512/';                     %updated 5/15/15
% % cd (fileDirectory);
% maxNumTrials = 150;
% tHold = 13;                                                                   %This is the hard-coded position that is reached to denote incorrect trials
% 
% % fileDirectory = '/Users/stubblefielde/Desktop/mice/mice2/';                 %for visitor re-approval trajectories
% 
%     [trialTimeVec, trialVec, numTrials, maxTrials, mousenumber, session] = psychFitRxnTime3LumsPos3(fileDirectory, lumVal_L1, lumVal_L2, lumVal_L3,...
%     lumVal_L4, lumVal_R1, lumVal_R2, lumVal_R3, lumVal_R4, tHold, maxNumTrials);
% 
% %% WORK ON THIS
%     trialNum = 1:(numTrials);    
%     newTimeMat = NaN(maxTrials, 1:trialNum);                                                
%     for i = 1:length(trialNum)
%         newTimeMat(1:maxTrials, i) = trialTimeVec(1:maxTrials, i) - trialTimeVec(1, i);
%     end                                                                         %These are new times aligned to trialStartTimes
%  
%     
% % Sanity check for trialVec = exact initial bar position :
%     initPos = trialVec(1,:)
% %     for p = 1:length(initPos)
% %         nineties = find(initPos(p)) == 90
% %     end
%     oneO5s = find(trialVec(1,:)) == abs(105)
%     oneTwenties = find(trialVec(1,:)) == abs(120)
%     oneThirty5s = find(trialVec(1,:)) == abs(135)
% %     end
% 
% %% Now plot position of bar with respect to time
% % For x-axis = time relative to the trial start time
% %     figure;
% %     hold on;
%     plot(newTimeMat(:,trialsToPlot1), trialVec(:,trialsToPlot1), 'k', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot2), trialVec(:,trialsToPlot2), 'bl', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot3), trialVec(:,trialsToPlot3), 'r', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot4), trialVec(:,trialsToPlot4), 'm', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot5), trialVec(:,trialsToPlot5), 'g', 'LineWidth',2);
%     plot(newTimeMat(:,trialsToPlot6), trialVec(:,trialsToPlot6), 'bl', 'LineWidth',2);
% %     plot(newTimeMat(:,trialsToPlot7), trialVec(:,trialsToPlot7), 'r', 'LineWidth',2);
% 
% 
%     
% %      plot(trialVec(:,trialsToPlot1), newTimeMat(:,trialsToPlot1), 'k', 'LineWidth',2);
% %      plot(trialVec(:,trialsToPlot2), newTimeMat(:,trialsToPlot2), 'g', 'LineWidth',2);
% % 
% %      plot(trialVec(:,trialsToPlot3), newTimeMat(:,trialsToPlot3), 'bl', 'LineWidth',2); %keep 
% %      plot(trialVec(:,trialsToPlot4), newTimeMat(:,trialsToPlot4), 'k', 'LineWidth',2);
% % 
% %      plot(trialVec(:,trialsToPlot5), newTimeMat(:,trialsToPlot6),'r', 'LineWidth',2);
% %      plot(trialVec(:,trialsToPlot6), newTimeMat(:,trialsToPlot5), 'm', 'LineWidth',2);
% 
% %     title(strcat(mousenumber,  session, '_ Unthresholded positions'));
% %     legend('trial 33','trial 50', 'trial 59')   %, 'trial 40', 'trial 50', 'location', 'NorthEast');
% %     set(gca, 'XLim', [0 1500]);
% %     set(gca, 'XTickLabel', [0:250:1500]);
% %     ylabel('Time from trial start (ms)', 'FontSize', 24);
% %     xlabel('(L)        Bar position         (R)', 'FontSize', 24);
% %     set(gca, 'XLim', [-300 300], 'FontSize', 20);
% %     set(gca, 'Ylim', [-150 150], 'FontSize', 20);
% %     hold off
% % 
% return