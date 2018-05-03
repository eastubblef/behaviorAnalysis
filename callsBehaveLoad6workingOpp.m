
function [taskbase] = callsBehaveLoad6workingOpp(fileDirectory, filenamestr)

fileDirectory = '/Users/stubblefielde/Desktop/behaviorFebMarAprMay16/Vgatseven/forAnalysis';
filenamestr = 'Vgatseven';

% if exist('filenamestr', 'var');
%     [filenamestrE, path] = uigetfile('*.csv*','select the csv file', fileDirectory);                 % get the actual # total trials 
% end
% if exist('filenamestr', 'var');
%     [filenamestrX, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fileDirectory);          % pxy file has solenoid discharge TS with respect to trial starts
% end

numfiles = 50;
rxPerformMat = []

dinfo = dir('*.csv*', '*pXY.csv');
for k = 1 : length(dinfo)
  thismat = dinfo(k).name;
  st = load(thismat);



   [taskbase] = behaveLoad6workingOpp(fileDirectory, filenamestr, csvFile, pXYfile)

end