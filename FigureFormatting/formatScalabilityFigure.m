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

methods = {'lag_under', 'lag_over', 'ccc', 'genzps'};

hf = figure();
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 1.5;

markstyle = {'^', 'x', 'o', 's'};
colors = {[0, 0, 1], [0, 0, 0], [1, 0, 0], [0, 1, 0]};
for lv = 1:length(methods)
    ndims = length(cmptimes.(methods{lv}).comptimes);
    hp = plot(2:ndims+1, cmptimes.(methods{lv}).comptimes, ...
        'Marker', markstyle{lv}, 'MarkerSize', 5, ...
        'Color', colors{lv}, 'MarkerFaceColor', colors{lv}, ...
        'MarkerEdgeColor', colors{lv});
    
    if strcmp(hp.Marker, 'x')
        hp.MarkerSize = 6;
    end
    
    % hold on / off
    if lv == 1
        hold on;
    elseif lv == length(methods)
        hold off;
    end
end

ha = gca;
ha.YScale = 'log';
ha.YLim = [10e-3, 10e4];
grid on;
ha.TickLabelInterpreter = 'latex';
ha.FontSize = FONT_SIZE;
ha.XLabel.String = 'Dimension';
ha.XLabel.FontSize = FONT_SIZE;
ha.XLabel.Interpreter = 'latex';
ha.YLabel.String = 'Computation Time [s]';
ha.YLabel.Interpreter = 'latex';
ha.YLabel.FontSize = FONT_SIZE;

lh = legend('\texttt{lag-under}', '\texttt{lag-over}', ...
    '\texttt{chance-open}', '\texttt{genzps-open}');
lh.Interpreter = 'latex';
lh.Position = [0.6332, 0.2159, 0.3206, 0.3582];

print('-dpng', '-r350', '../exampleFigs/pngs/scalability-example.png');
