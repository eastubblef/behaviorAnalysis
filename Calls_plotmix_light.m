
function [percent_left_odor1, fit_values1, mix_info1, percent_left_odor2, fit_values2, mix_info2...
   percent_left_odor3, fit_values3, mix_info3] = Calls_plotmix_light(rbase_fname1, rbase_fname2, rbase_fname3);


% This wrapper function calls PlotMixPerformance.m in order to evaluate
% behavioral sessions in which drug was delivered and compare to sessions in which saline was delivered


%% Input

rbase_fname1 = 'C:\Felsen Lab docs\Matlab\Data\behavior\Plots\ChatCreChR1\tb_ChatCreChR1_131105a'
% saline_1 = 'tb_WT1_131007a_saline'

%For calling up all trials:
load(rbase_fname1);
num_trials = numel(taskbase.stimID);
light_on_trials = find(taskbase.fib_lit == 1);
light_off_trials  = find(taskbase.fib_lit == 0);


%% Call it up for light_off trials:

fn = rbase_fname1

for trials_to_skip = 0
    mpmp_flag = 0;
    divi = 0;
    no_plot_flag = 1;
    
    
    [percent_left_odor1, fit_values1, mix_info1] = PlotMixPerformance(rbase_fname1, no_plot_flag, trials_to_skip, mpmp_flag, divi);
    
    % Plot PlotMixPerformance graph as normal
        tot_trials = mix_info1.num_Lf + mix_info1.num_Rf;
        num_Lf = mix_info1.num_Lf 
        Lf_trials = (num_Lf) ./ (tot_trials);
%         saline_pre = percent_left_odor1; %To help set up legend at end

        p = plot(saline_pre, Lf_trials, 'k.');
        hold on;

        set(p, 'MarkerFaceColor', 'k', 'MarkerSize', 20);
        set(gca, 'YLim', [0 1], 'XLim', [0 100]);

        % the fit
        p = plot(fit_values1.x_axis, fit_values1.y_fit, 'k');
        set(p, 'LineWidth', 2);

        % Labeling commands
        xlabel('% Odor A', 'FontSize', 16);
        ylabel('Fraction Left choices', 'FontSize', 16);
        title([num2str(rbase_fname1) '_saline']);


        plot(fit_values1.x_axis, fit_values1.y_fit,  'k-',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerSize', 15);
%     plot(percent_left_odor1, mix_info1.num_Lf, 'k.', 'MarkerSize', 20);
%     hold on;
%     xlabel('Percent Odor A'); ylabel('Fraction Left Choices');
%     hold on;
end

%% Drug plot (red) to overlay:

rbase_fname2 = 'C:\Felsen Lab docs\Matlab\Data\behavior\Plots\ChatCreChR1\tb_ChatCreChR1_131103a'
string2 = '[]'
fn = rbase_fname2;
   
for trials_to_skip = 0;
    mpmp_flag = 0;
    divi = 1;
    no_plot_flag = 1;
    
    [percent_left_odor2, fit_values2, mix_info2] = PlotMixPerformance(rbase_fname2, no_plot_flag);
    
     hold on;
    % Plot PlotMixPerformance graph as normal for drug day
            
        tot_trials2 = mix_info2.num_Lf + mix_info2.num_Rf;
        num_Lf2 = mix_info2.num_Lf 
        Lf_trials2 = (num_Lf2) ./ (tot_trials2);
        DHBE = percent_left_odor2;
        p2 = plot(DHBE, Lf_trials2, 'r.')

        hold on;
        set(p2, 'MarkerFaceColor', 'k', 'MarkerSize', 10);
        set(gca, 'YLim', [0 1], 'XLim', [0 100]);

        % the fit
        p2 = plot(fit_values2.x_axis, fit_values2.y_fit, 'r');
        set(p2, 'LineWidth', 2);

        % Labeling commands
        xlabel('% Odor A', 'FontSize', 16);
        ylabel('Fraction Left choices', 'FontSize', 16);
        title([num2str(rbase_fname2) '_DHBE']);


%     plot(fit_values2.x_axis, fit_values2.y_fit,  'r-',...
%                 'LineWidth',2,...
%                 'MarkerEdgeColor','r',...
%                 'MarkerSize', 15);
%     
      %plot(percent_left_odor2, mix_info2.num_Lf, 'r.', 'MarkerSize', 20);
    
      hold on;
end
    
%% For 2nd saline session "Post-saline":

rbase_fname3 = 'C:\Felsen Lab docs\Matlab\Data\behavior\Plots\ChatCreChR1\tb_ChatCreChR1_131102a'

fn3 = rbase_fname3
 
for trials_to_skip = 0
    mpmp_flag = 0;
    divi = 0;
    no_plot_flag = 1;

    
    [percent_left_odor3, fit_values3, mix_info3] = PlotMixPerformance(rbase_fname3, no_plot_flag);
      
    % Plot PlotMixPerformance graph as normal
        tot_trials = mix_info3.num_Lf + mix_info3.num_Rf;
        num_Lf = mix_info3.num_Lf 
        Lf_trials = (num_Lf) ./ (tot_trials);
        saline_post = percent_left_odor3; %Trying to get legend to work

        p = plot(saline_post, Lf_trials, 'b*');

        hold on;

        set(p, 'MarkerFaceColor', 'k', 'MarkerSize', 10);
        set(gca, 'YLim', [0 1], 'XLim', [0 100]);

        % the fit
        p = plot(fit_values3.x_axis, fit_values3.y_fit, 'k-');
        set(p, 'LineWidth', 2);

        % Labeling commands
        xlabel('% Odor A', 'FontSize', 16);
        ylabel('Fraction Left choices', 'FontSize', 16);
        title([num2str(rbase_fname3) '_saline']);
        plot(fit_values3.x_axis, fit_values3.y_fit,  'b-',...
            'LineWidth',2,...
            'MarkerEdgeColor','g',...
            'MarkerSize', 15);
        hL = legend ('salinepre', 'DHBE red', 'salinepost blue', 'location', 'SouthEast');

end


end