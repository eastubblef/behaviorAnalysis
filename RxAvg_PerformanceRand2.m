
%% This script is for generating a heat map (imagesc) of binned rx times' %corrects 
% Bin sizes are set for 200 ms = fastest bin; 200-600 ms = fast bin; 600-1000 = mid bin; >1000 = slow bin
% 11.2017  This script calls rxSort for forming the matrix
%% Inputs:
fpath = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
cd(fpath);
% folder = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/4image_tbRand';
% fullfile = strcat(folder,'/',file);
numfiles = 70;
rxPerformMat = []

dinfo = dir('*.mat');
for k = 1 : length(dinfo)
  thismat = dinfo(k).name;
  st = load(thismat);
  rxMat = st.taskbase.rxMat;
  
   [rxSort] = rxSort4Mat(rxMat);
   rxSort = rxSort';
   
%    [Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000] = rxSortMat4perform2(rxSort);
   [Num100, Num200, Num400, Num600, Num1000, Num1000to3000] = rxSortMat4perform2(rxSort);

%     NumFastest = NumFastest;
%     NumFast = NumFast;
%     NumMid = NumMid;
%     NumSlow = NumSlow;
%     cat = horzcat(Num100, Num200, Num300, Num400, Num500, Num600, Num700, Num800, Num900, Num1000, Num1100, Num1200, Num1300, Num1400, Num1500, Num1600, Num1700, Num1800, Num1900, Num2000);
    cat = horzcat(Num100, Num200, Num400, Num600, Num1000, Num1000to3000);

    rxPerformMat(k,:) = cat;

   clear rxSort;
   clear cat;
end

%% Create bins for rx times and insert % corrects
% catPerform = vertcat(cat1,cat2,cat3,cat4);

h_image = imagesc(rxPerformMat);
colorbar;
ylabel('Randomized sessions'); xlabel('Rx times *100 (ms)'); 

