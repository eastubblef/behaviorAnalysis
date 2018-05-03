%% AnalysesFigSetUp.m
%% Sets default values for figures.
%% For general help on object properties, in helpdesk search
%% "Function Name" for "get", and at bottom of page,
%% click on "Handle Graphics Properties"
%%
%% For more info on defaults, in helpdesk search "Document Titles"
%% for "Defining Default Values" and "Setting Default Property Values".
%% 10/27/03
%%
%% to re-reference entire figure, after calling this fcn: set(gcf, 'CurrentAxes', whole_fig);
%%
%% $ Settings updated for selectivity plots for rat neurons 2/13/06 - GF $

GetSCGlobals;
global RESOLUTION; % set to 1000 to plot spike times to 1 msec resolution
global TRIAL_EVENTS;

%PRESENTATION_TYPE = 'talk'; % talk or paper (affects font sizes)
PRESENTATION_TYPE = 'paper'; % talk or paper (affects font sizes)

%% set a square figure window, so that when units are normalized, .1 vertically = .1 horizontally

% get width of current monitor
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);

WIDTH_FRACTION = 0.555555555555556; % based on what worked on laptop

FIG_WIDTH = screen_width * WIDTH_FRACTION;
%FIG_WIDTH = 660; % good for desktop (Slothrop)
%FIG_WIDTH = 800; % good for laptop

FIG_HEIGHT = FIG_WIDTH;
fig_position = [350 20 FIG_WIDTH FIG_HEIGHT];

PANEL_LABEL_FONT_SIZE = 20; %% size of the letters labeling each subfigure panel
PANEL_LABEL_FONT_WEIGHT = 'normal';
%PANEL_LABEL_FONT_WEIGHT = 'bold';

%AXIS_LABEL_FONT_SIZE = 13; %% size of axis labels
%TEXT_FONT_SIZE = 13; %% size of axis labels

if strcmp(PRESENTATION_TYPE, 'talk')
    TEXT_FONT_SIZE = 16; %% size of axis labels
else
    %TEXT_FONT_SIZE = 12; %% size of axis labels
    
    %TEXT_FONT_SIZE = 14; %% size of axis labels - corresponds to 11.97 in AI CS3
    TEXT_FONT_SIZE = 12; %% changed 12/2/09 for Slothrop (office desktop)

end
CAP_OFFSET = 0.005; %% "error" between top of .AI objects and vertical alignment of text with "cap" option
DIST = 80; % cm from cat to monitor - necessary for sizing scale bars

fig = figure;
set(fig, 'Position', fig_position);

set(fig,...
    'Units', 'normalized',...
    'DefaultAxesFontName', 'Arial',...
    'DefaultAxesColor', [1 1 1],...
    'DefaultAxesFontSize', TEXT_FONT_SIZE,...
    'DefaultAxesLineWidth', 1,...
    'DefaultAxesColorOrder', [0 0 0],...
    'DefaultAxesTickDir', 'out',...
    'DefaultAxesTickLength', [0.015 0.015],...
    'DefaultLineLineWidth', 2,...
    'DefaultTextFontName', 'Arial',...
    'DefaultTextFontSize', TEXT_FONT_SIZE, ...
    'DefaultTextVerticalAlignment', 'Middle', ...
    'DefaultTextHorizontalAlignment', 'Center' ...
    );

% set a 'whole_fig' axis that can be rereferenced in the figure-producing
% mfile to make axis labels, etc.
whole_fig = axes('Position', [0 0 1 1]);

if exist('SHOW_GRIDS')
    if SHOW_GRIDS == 1
        set(gca,'Visible','on');
        grid on;
    else
        set(gca,'Visible','off');
    end
else
    set(gca,'Visible','off');
end    

%%% Below is specific for making rasters/psths from ratbase files %%%

%% Define the trial events of interest - correspond to fields of ratbase structure.
%% Get these from the GetSCGlobals definitions, but rename some to be
%% compatible with ratbase fieldnames.
%% Note that these are not chronological, b/c GoToneOn was added later.
trial_events = fieldnames(TRIAL_EVENTS);
trial_events{2} = 'DIO'; %DIO is OdorValveOn
trial_events{5} = 'WaterDeliv'; %WaterDeliv is WaterValveOn


%% Set the colors common to rasters and psths
RASTER_COLORS.L = [0 0 0];
RASTER_COLORS.R = [1 0 0];
RASTER_COLORS.mix.odorA = [0 0 1; 0.5 0.5 0.5];
RASTER_COLORS.mix.odorB = [0 1 0; 0.5 0 0.5];

RASTER_COLORS.correct = '-';
RASTER_COLORS.error = '--';

RASTER_COLORS.gotone_short = '-';
RASTER_COLORS.gotone_long = '--';

% colors of raster ticks for trial events
SECONDARY_EVENTS_COLORS = [...
    0 1 0;...
    1 0.5 1;...
    0.75 0 0.75;...
    0 0.75 0.75;...
    0 0 1;...
    0.75 0 0;...
    0 1 0;...
    0.85 0.85 0];



