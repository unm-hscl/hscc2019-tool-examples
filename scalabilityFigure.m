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
cmptimes = load('scalability_comptimes.mat');

methods = {'lag', 'ccc', 'genzps'};

hf = figure();
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 1.5;

for lv = 1:length(methods)
    ndims = length(cmptimes.(methods{lv}).comptimes);
    semilogy(2:ndims+1, cmptimes.(methods{lv}).comptimes);
    
    % hold on / off
    if lv == 1
        hold on;
    elseif lv == length(methods)
        hold off;
    end
end

ha = gca;
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

lh = legend('\texttt{lag-under}', '\texttt{chance-open}', ...
    '\texttt{genzps-open}');
lh.Interpreter = 'latex';
