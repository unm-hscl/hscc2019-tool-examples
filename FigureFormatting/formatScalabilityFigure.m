%
% Name: scalabilityFigure.m
% Author: Joseph D. Gleason
%
% Description: Plot the scalability figure for the HSCC paper
%

close all;
clearvars;

% paper column width in pts
HSCC_COLWIDTH = 241.14749;

FONT_SIZE = 8;

% load the data
cmptimes = load('../MatFiles/scalability_comptimes.mat');

methods = {'ccc', 'genzps'};

hf = figure();
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 1.5;

markstyle = {'^', 's', 'o'};
colors = {[0, 0, 1], [0, 0, 0], [1, 0, 0]};
for lv = 1:length(methods)
    ndims = length(cmptimes.(methods{lv}).comptimes);
    hp = plot(2:ndims+1, cmptimes.(methods{lv}).comptimes, ...
        'Marker', markstyle{lv}, 'MarkerSize', 5, ...
        'Color', colors{lv}, 'MarkerFaceColor', colors{lv}, ...
        'MarkerEdgeColor', colors{lv});
        
    % hold on / off
    if lv == 1
        hold on;
    elseif lv == length(methods)
        hold off;
    end
end
hold on;
load('../../../2017/CSSL/Figure1_curseOfDim_full_10_5_20points_n1to40_reqdOnly.mat', 'elapsed_time_DP')
plot(2:4,[elapsed_time_DP{:}],'ro-','MarkerFaceColor', 'r','DisplayName','Dynamic prog.');

ha = gca;
ha.YScale = 'log';
ha.YLim = [10^(-1), 10^(4.5)];
grid on;
ha.TickLabelInterpreter = 'latex';
ha.FontSize = FONT_SIZE;
ha.XLabel.String = 'Dimension';
ha.XLabel.FontSize = FONT_SIZE;
ha.XLabel.Interpreter = 'latex';
ha.YLabel.String = 'Computation Time (s)';
ha.YLabel.Interpreter = 'latex';
ha.YLabel.FontSize = FONT_SIZE;

lh = legend('\texttt{chance-open}', ...
    '\texttt{genzps-open}','Dynamic programming','Location','SouthEast');
lh.Interpreter = 'latex';
