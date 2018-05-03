%% PlotMixPerformance
% Plots performance on mixture discrimination, for 1 mixture pair.
%
% USAGE: [percent_left_odor, percent_L, fits, mix_info] = PlotMixPerformance(rbase_fname, no_plot_flag, trials_to_skip, mpmp_flag, divi)
% EXAMPLE: PlotMixPerformance
%
% INPUTS:  rbase - the ratbase filename (created by solo2ratbase) w/ path, if not in current folder
%              default - user selects ratbase file w/ gui
%          no_plot_flag - when nonzero, data is not plotted
%          trials_to_skip - number of trials at beginning of session to
%              disregard (b/c rat hasn't "gotten into it" yet) (OPTIONAL,
%              but note: if using this variable, all previous inputs must
%              also be entered)
%          mpmp_flag - added to allow the subplot function run properly in
%          the MultiPlotMixPerformance (mpmp) mfile
%          divi - '0' means the user does not wish to divide the trial
%          information in taskbase into two portions. '1' means they do
%          wish to do this (allows comparison with fiber light active
%          'Mix2afc' information)
%
% OUTPUTS: percent_left_odor - 1 x num_mixes vector of percent of odor A
%          percent_L - 1 x num_mixes vector of percent left choice
%          pair_num - 1 x num_mixes vector of which odor pair each stim is
%              calculated bias
%          mix_info - other info to return, such as:
%              num_L - total number of L choices for each mixture ratio
%              num_R - total number of R choices for each mixture ratio
% 
% Note that there may be 0 trials for a particular mix_ind (if so,
% percent_L will be NaN).
%
% 3/19/07 - GF
% $ Update 8/1/07 GF - variables are outputted
% $ Update 4/30/08 GF - trials_to_skip optional input added
% $ Update 08/21/2011 JDC - Modified previous mfile to accept the
%    'taskbase' variable.  Also, this version removes the ability to plot more
%    than 1 odor mixture pair
% $ Updated 10/20/2011 JDC - Added code to allow signmodial plots for the
%     optogenetic experiments to be generated (i.e. 'fib_lit' variable in
%     taskbase).  Separated the code that creates the plot into a
%     subfunction called 'sigmo_plot'
% $ Updated 11/25/2011 JDC - Added code to randomly split trials into two
%   plots (using the input variable 'divi') in order to compare the
%   fiber light experiment plots with regular task days.

% GF (8/26/13) noticed 'percent_L' not assigned:
% function [percent_left_odor, percent_L, fit_values, mix_info] = PlotMixPerformance(rbase_fname, no_plot_flag, trials_to_skip, mpmp_flag, divi)

function [percent_left_odor, fit_values, mix_info] = PlotMixPerformance(rbase_fname, no_plot_flag, trials_to_skip, mpmp_flag, divi)

if nargin < 1 || isnumeric(rbase_fname)
    
    if nargin == 1
        no_plot_flag = rbase_fname;
    else
        no_plot_flag = 0;
    end
    
    [fn pn] = uigetfile('C:\behavior\ExperPort\Data\Data\Jamie\*.mat','Select a file to transform');
    rbase_fname = [pn fn];
    
    % If the uigetfile screen is escaped from...
    if fn == 0,
        percent_left_odor = '';
        percent_L = '';
        fit_values = '';
        mix_info = '';
        return;
    end;

elseif nargin < 2 % the taskbase filename was entered, but not a no_plot_flag
    
    no_plot_flag = 0;
    
    fn = rbase_fname; % need to know the filename
else 
    fn = rbase_fname; % need to know the filename
end


load(rbase_fname);

%% remove the first trials_to_skip trials
if nargin > 2,% trials_to_skip was entered (otherwise, don't skip any trials)    
    fnames = fieldnames(taskbase);
    for fname_ind = 1:length(fnames)
        if isnumeric(taskbase.(fnames{fname_ind}))            
            taskbase.(fnames{fname_ind}) = taskbase.(fnames{fname_ind})((trials_to_skip + 1):end);
        end
    end
end

%% Determines if the 'mpmp_flag' and 'divi' variables are present, if not it sets the value
if ~ismember(who, 'mpmp_flag'), 
    mpmp_flag = 0;
end

if ~ismember(who, 'divi'), 
    divi = 0;
end

%% This portion of code extracts the percent left values for each stimulus

% Calculates the number of odorants used
mix_ind = find(~isnan(cell2mat(taskbase.odor_params.odors)));
num_mixes = length(mix_ind);  

% The prescence of 'n' or 'f' was added because of the 'Mix2afc_Light' (or
% 'M2L') protocol. 'n' -> fiber light on, 'f' -> fiber light off.  For the 
% 'Mix2afc' protocol these are simply ignored.  10/09/2011.
    
% This if-else structure calculates the fraction of left choice and right
% choices for a 'Mix2afc' task.  The main if-else structure is concerned
% with the presence of a fiber light in the trial information.
if ~isfield(taskbase, 'fib_lit') | isnan(taskbase.fib_lit) % Means there is no fiber light present in trial information
    if divi == 0,  % Random division in trial information not wanted
        num_Ln = NaN;
        num_Rn = NaN;
        num_Lf = 99 * ones(1, num_mixes);
        num_Rf = 99 * ones(1, num_mixes);

        % Determiens %choice L and R for each stimulus
        for k = 1:num_mixes    
            num_Lf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 1)));   
            num_Rf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 2))); 
        end

        %percent_L = (num_L ./ (num_L + num_R)) * 100;
        percent_Ln = NaN;
        percent_Lf = (num_Lf ./ (num_Lf + num_Rf));
    
    elseif divi == 1, % User wishes to randomly separate trial
        num_Ln = 99 * ones(1, num_mixes);
        num_Rn = 99 * ones(1, num_mixes);
        num_Lf = 99 * ones(1, num_mixes);
        num_Rf = 99 * ones(1, num_mixes);
        % Variable to separate trials into two groups
        rando = rand(1, numel(taskbase.choice));

        % Determiens %choice L and R for each stimulus
        for k = 1:num_mixes    
            num_Ln(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 1) & (rando' >= 0.5)));   
            num_Rn(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 2) & (rando' >= 0.5)));
            num_Lf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 1) & (rando' < 0.5)));   
            num_Rf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 2) & (rando' < 0.5)));    
        end

        %percent_L = (num_L ./ (num_L + num_R)) * 100;
        percent_Ln = (num_Ln ./ (num_Ln + num_Rn));
        percent_Lf = (num_Lf ./ (num_Lf + num_Rf));    
    end
    
else % Fiber light present, and used
    num_Ln = 99 * ones(1, num_mixes);
    num_Rn = 99 * ones(1, num_mixes);
    num_Lf = 99 * ones(1, num_mixes);
    num_Rf = 99 * ones(1, num_mixes);
    
    % Determiens %choice L and R for each stimulus
    for k = 1:num_mixes    
        num_Ln(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 1) & (taskbase.fib_lit' == 1)));   
        num_Rn(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 2) & (taskbase.fib_lit' == 1)));
        num_Lf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 1) & (taskbase.fib_lit' == 0)));   
        num_Rf(k) = length(find((taskbase.stimID == mix_ind(k)) & (taskbase.choice == 2) & (taskbase.fib_lit' == 0)));    
    end
    
    %percent_L = (num_L ./ (num_L + num_R)) * 100;
    percent_Ln = (num_Ln ./ (num_Ln + num_Rn));
    percent_Lf = (num_Lf ./ (num_Lf + num_Rf));
end

%% Sets the odorant concentration variable 'percent_left_odor'

% Pulls odorant concentration from 'taskbase' variable
percent_left_odor = cell2mat(taskbase.odor_params.odors(mix_ind));

% Alters the odorant concentration variable ('percent_left_odor') from the 
% to reflect the increasing concentrations with respect to left choices
for k = 1:length(percent_left_odor),
    if percent_left_odor(k) > 0,
        percent_left_odor(k) = 100 - percent_left_odor(k);
    end
end

% Eliminates "negative" concentrations ("negative" concentrations simply
% denoted the chirality of the enantiomer [i.e. +/- carvone])
percent_left_odor = abs(percent_left_odor);


%% get the best binomial fit, but only w/ the mixture pair! (bug fix 3/21/08)
% note that mix ratios w/ 0 trials (i.e., the default 95s) do not
% contribute to the fit, so it's ok that they're included in the
% percent_left_odor

% 'f' is for light off, 'n' is for light on
if isnan(num_Ln),
    bf = glmfit(percent_left_odor', [num_Lf' (num_Lf + num_Rf)'], 'binomial');

    x_axisf = 0:100;
    y_fitf = glmval(bf, x_axisf, 'logit');

    % return the fit y and the x it was based on, and the calculated bias
    fit_values.x_axis = x_axisf;
    fit_values.y_fit = y_fitf;
    fit_values.bias = 50 + (bf(1)/bf(2)); % Bobby's calculation, based on Hatim's email explanation
else
    bf = glmfit(percent_left_odor', [num_Lf' (num_Lf + num_Rf)'], 'binomial');

    x_axisf = 0:100;
    y_fitf = glmval(bf, x_axisf, 'logit');

    % return the fit y and the x it was based on, and the calculated bias
    fit_values.x_axisf = x_axisf;
    fit_values.y_fitf = y_fitf;
    fit_values.biasf = 50 + (bf(1)/bf(2)); % Bobby's calculation, based on Hatim's email explanation
    
    bn = glmfit(percent_left_odor', [num_Ln' (num_Ln + num_Rn)'], 'binomial');

    x_axisn = [0:100];
    y_fitn = glmval(bn, x_axisn, 'logit');

    % return the fit y and the x it was based on, and the calculated bias
    fit_values.x_axisn = x_axisn;
    fit_values.y_fitn = y_fitn;
    fit_values.biasn = 50 + (bn(1)/bn(2)); % Bobby's calculation, based on Hatim's email explanation    
end;


%% package additional info structure
mix_info.num_Ln = num_Ln;
mix_info.num_Rn = num_Rn;
mix_info.num_Lf = num_Lf;
mix_info.num_Rf = num_Rf;

%% Plot function execution

% Sends info to subfunction to create a more clear script

light = 0;
sigmo_plot(fn, no_plot_flag, percent_left_odor, percent_Lf, x_axisf, y_fitf, mpmp_flag, light)

% Executes for 'sigmo_plot' function if the protocol employed a fiber light
if ~isnan(num_Ln),
    fn = [fn(1:end-5), '-On', fn(end-4:end)];
    light = 1;
    sigmo_plot(fn, no_plot_flag, percent_left_odor, percent_Ln, x_axisn, y_fitn, mpmp_flag, light)
end;
        

end


function [] = sigmo_plot(fn, no_plot_flag, percent_left_odor, percent_L, x_axis, y_fit, mpmp_flag, light)
    %% Plotting commands

    if no_plot_flag == 0

    % The 'mpmp_flag' was introduce to allow plotting performance
    % across days using a subplot function (JDC 10/13/2011)
        if mpmp_flag ~= 0,
            hold on;

            % 1st pair
            if light
                p = plot(percent_left_odor, percent_L, 'r+');
                set(p, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
            else
                p = plot(percent_left_odor, percent_L, 'ko');
                set(p, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
            end

            % the fit
            if light
                p = plot(x_axis, y_fit, 'r');
                set(p, 'LineWidth', 1);
            else
                p = plot(x_axis, y_fit, 'k');
                set(p, 'LineWidth', 1);
            end;

            % Labeling commands - No labeling commands because of
            % congestion
            %xlabel('% Odor A', 'FontSize', 5);
            %ylabel('Fraction Left choices', 'FontSize', 5);
            
            set(gca, 'YLim', [0 1], 'XLim', [0 100]);
        else % 'If the MultiPlotMixPerformance' mfile has not been executed
    % Plots PlotMixPerformance graph as normal
            figure;
            hold on;

            % 1st pair
            p = plot(percent_left_odor, percent_L, 'ko');
            set(p, 'MarkerFaceColor', 'k', 'MarkerSize', 10);

            % the fit
            p = plot(x_axis, y_fit, 'k');
            set(p, 'LineWidth', 2);

            % Labeling commands
            xlabel('% Odor A', 'FontSize', 16);
            ylabel('Fraction Left choices', 'FontSize', 16);

    % Added Code to title the graph properly 
            % the following are where you define information for the new file names
            rat = fn((end-15):(end-5));              % rat (assume 4 characters)

            % A portion of code that necessitates the variable 'rat' has an '_'
            % in it.  This allows the entire name of the 'rat' to be determined
            k = 1;
            while strcmp(rat(1), '_') ~= 1 % actually, this is a 3 character-named rat
                rat = fn((end-(15 + k)) : (end-5));
                k = k + 1;
            end

            % Removes the underscore ('_')
            rat = rat(2:end); 
            % Removes underscore (because the funciton 'title' has difficulty
            % with it)
            rat(strfind(rat, '_')) = '-';

    % Titles graph        
            title(rat);
            set(gca, 'YLim', [0 1], 'XLim', [0 100]);
        end
    end
end