%% This script calls excel2taskbase_bar for behavioral extraction
%   BS 2.29.16
%   BS 7/3/15

%   INPUT
%   rawdatafile = .csv file 
%   th = threshold (for "correct" trial)
%   bound = maximum eccentricity of the bar 
%   
%   OUTPUT
%   taskbase_bar:

fpath = '/Users/stubblefielde/Desktop/mfiles/DudmanLab/data/spikes/recordings15/Vgatthree/2ndpass/better/behaveChunks'
% fpath = '/Users/stubblefielde/Desktop/VGATone/062315lums/furtherPositions'
% newPath = strcat(fpath, '/', 'Gad2_tag', '/', '150529ipsi', '/');

if exist('fpath', 'var');
%     [rawdatafile, path] = uigetfile('*.csv*','select the .csv file', fpath)            
      [rawdatafile, path] = uigetfile('*pXY.csv*','select the pXY.csv file', fpath);            
end

% 2.28.16 Now with 2AFC task, bar must be centered; else, incorrect trial -
% define in the core function since this is reported in the pxy file
% bound = 200;                                                                  % Maxiumum eccentricity of bar
% th = bound;                                                                   % For "correct" movement detection

% [taskbase_bar] = excel2taskbase_bar (rawdatafile , th , bound);
[taskbase_bar] = excel2taskbase_bar (rawdatafile);

initiation = taskbase_bar.initiation;                                       % Time (s) of 'choice' movement start

PresentationFigSetUp;
hold on;
subplot(2,2,1);
x = 1:numel(initiation);
initPlot = scatter(x, initiation);
ylabel('Initiation time (s)');
xlabel('Trial number');
% xlabels = [0 2 4 8 16]
set(gca, 'XTick', [0:10:length(initiation)], 'Ylim', [0 2.5]);
set(gca, 'XTickLabel', [0:10:numel(initiation)], 'Xlim', [0:numel(initiation)]);
