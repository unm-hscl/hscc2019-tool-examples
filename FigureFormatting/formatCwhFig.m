% 
% Name:        formatCwhFig.m
% Author:      Joseph D. Gleason
% 
% Description: Format the CWH example figure for printing
% 

close all;
clearvars;

% paper column width in pts
HSCC_COLWIDTH = 241.14749;

FONT_SIZE = 8;

% Load figures
figs(2) = openfig('../exampleFigs/CWH_example_ZeroInitVel_traj.fig');
figs(1) = openfig('../exampleFigs/CWH_example_ZeroInitVel.fig');

% remove the legend
delete(figs(1).Children(1));

% get axes
axs(2) = figs(2).Children(1);
axs(1) = figs(1).Children(1);

% delete / copy stuff from one figure to the other
delete(axs(1).Children(1:2));
delete(axs(2).Children(length(axs(2).Children)));
delete(axs(2).Children(length(axs(2).Children)));
delete(axs(2).Children(length(axs(2).Children)));

copyobj(axs(2).Children, axs(1));
delete(figs(2));

hf = figs(1);
ha = axs(1);

clearvars figs axs;

% undock the figures
hf.WindowStyle = 'normal';

% need to pause to let MATLAB catch up to the undocking
pause(2);

% reposition figures
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 3.2;
hf.Position(1) = 3;
hf.Position(2) = 3;

% get rid of ratioed scaling
ha.DataAspectRatioMode = 'auto';
ha.PlotBoxAspectRatioMode = 'auto';

ha.FontSize = FONT_SIZE;
ha.TickLabelInterpreter = 'latex';
ha.XLabel.Interpreter = 'latex';
ha.XLabel.FontSize = FONT_SIZE;
ha.YLabel.Interpreter = 'latex';
ha.YLabel.FontSize = FONT_SIZE;
ha.Units = 'inches';
ha.Position(4) = 2;
ha.Units = 'normalized';
ha.Position(2) = 0.34;
ha.YLim = [-1.5, 0.1];
ha.XLim = [-1.5, 1.5];

delete(ha.Children(1));
delete(ha.Children(7));

for lv = 1:length(ha.Children)
    if strcmp(ha.Children(lv).Type, 'scatter')
        ha.Children(lv).LineWidth = 1;
        ha.Children(lv).SizeData = 15;
    elseif strcmp(ha.Children(lv).Type, 'line')
        ha.Children(lv).MarkerSize = 3;
        ha.Children(lv).LineWidth = 1;
        ha.Children(lv).Color = [0, 0.8, 0];
        ha.Children(lv).MarkerFaceColor = [0, 1, 0];
        ha.Children(lv).MarkerEdgeColor = [0, 1, 0];
    end
    
    switch(lower(ha.Children(lv).DisplayName))
        case 'lagrangian'
            ha.Children(lv).FaceColor = [0, 0.9, 0];
            ha.Children(lv).FaceAlpha = 1;
        case 'fourier transform'
            ha.Children(lv).FaceColor = [0, 0.6, 1];
        case 'target set'
            ha.Children(lv).FaceColor = [0, 0, 0];
        case 'safe set'
            ha.Children(lv).FaceColor = [0.95, 0.95, 0];

            ha.Children(lv).Vertices = [1.5, -1.5, 0; ...
                0, -0, 0; ...
               -1.5, -1.5, 0];
        case 'chance constraint'
            ha.Children(lv).FaceColor = [1, 0.6, 0];
        case 'bad trajectory'
            ha.Children(lv).Color = [1, 0, 0];
            ha.Children(lv).LineStyle = '-';
            ha.Children(lv).LineWidth = 1.5;
            ha.Children(lv).Marker = 'x';
            ha.Children(lv).MarkerSize = 7;
            ha.Children(lv).MarkerFaceColor = [1, 0, 0];
            ha.Children(lv).MarkerEdgeColor = [1, 0, 0];
        otherwise
            % do nothing
    end
end


% annotations galore!!!
% Safe set
% -----------------------
% the box
% an = annotation(hf, 'textbox');
% an.String = 'Sets';
% an.Interpreter = 'latex';
% an.FontSize = FONT_SIZE;
% an.Position(1) = 0.05;
% an.Position(2) = ha.Position(2) - 0.15;
% an.Position(4) = 0.05;
% an.Position(3) = 0.3;
% an.LineStyle = 'none';
% pos = an.Position;

delete(ha.Children(9));

an = annotation(hf, 'rectangle');
an.Position(1) = 0.1;
an.Position(2) = ha.Position(2) - 0.18;
an.Units = 'inches';
an.Position(3) = 0.1;
an.Position(4) = 0.1;
an.Units = 'normalized';
an.FaceColor = [0.95, 0.95, 0];

% Text
txt = annotation(hf, 'textbox');
txt.Position(1) = an.Position(1) + 0.04;
txt.Position(2) = an.Position(2) - 0.002;
txt.Position(3) = 0.3;
txt.Position(4) = 0.05;
txt.String = 'Safe Set';
txt.Interpreter = 'latex';
txt.LineStyle = 'none';
txt.FontSize = FONT_SIZE;

strs = {'target', 'ccc', 'ftgenz'};
set_colors = {[0, 0, 0], [1, 0.6, 0], [0, 0.6, 1]};
lefttab = an.Position(1);
righttab = 0.5;
for lv = 1:length(strs)
    if mod(lv, 2) == 1
        % odd
        % the box
        pos = an.Position;
        an = annotation(hf, 'rectangle');
        an.Position(1) = righttab;
        an.Position(2) = pos(2);
        an.Position(3) = pos(3);
        an.Position(4) = pos(4);
        an.Units = 'normalized';
        an.FaceColor = set_colors{lv};
        if strcmp(strs{lv}, 'ccc')
            an.FaceColor = set_colors{lv-1};
            an = annotation(hf, 'rectangle');
            an.Position(1) = righttab;
            an.Position(2) = pos(2);
            an.Position(3) = pos(3);
            an.Position(4) = pos(4);
            an.Units = 'normalized';
            an.FaceColor = set_colors{lv};
            an.FaceAlpha = 0.9;
        end
    else
        % even
        % the box
        pos = an.Position;
        an = annotation(hf, 'rectangle');
        an.Position(1) = lefttab;
        an.Position(2) = pos(2) - 0.05;
        an.Position(3) = pos(3);
        an.Position(4) = pos(4);
        an.Units = 'normalized';
        an.FaceColor = set_colors{lv};
        if strcmp(strs{lv}, 'ccc')
            an.FaceColor = set_colors{lv-1};
            an = annotation(hf, 'rectangle');
            an.Position(1) = lefttab;
            an.Position(2) = pos(2) - 0.05;
            an.Position(3) = pos(3);
            an.Position(4) = pos(4);
            an.Units = 'normalized';
            an.FaceColor = set_colors{lv};
            an.FaceAlpha = 0.9;
        end
    end
   

    % Text
    txt = annotation(hf, 'textbox');
    txt.Position(1) = an.Position(1) + 0.04;
    txt.Position(2) = an.Position(2) - 0.002;
    txt.Position(3) = 0.3;
    txt.Position(4) = 0.05;
    switch(strs{lv})
        case 'lagrangian'
            txt.String = '\texttt{lag-under}';
        case 'dynprog'
            txt.String = '\texttt{SReachDynProg}';
        case 'ccc'
            txt.String = '\texttt{chance-open}';
        case 'ftgenz'
            txt.String = '\texttt{genzps-open}';
        case 'target'
            txt.String = 'Target Set';
        otherwise
            % do nothing
    end 
    txt.Interpreter = 'latex';
    txt.LineStyle = 'none';
    txt.FontSize = FONT_SIZE;
end

% pos = an.Position;

% an = annotation(hf, 'textbox');
% an.String = 'Trajectories';
% an.Interpreter = 'latex';
% an.FontSize = FONT_SIZE;
% an.Position(1) = 0.05;
% an.Position(2) = pos(2) - 0.08;
% an.Position(4) = 0.05;
% an.Position(3) = 0.3;
% an.LineStyle = 'none';
% pos = an.Position;

% ax = axes(hf);
% ax.Position(1) = 0.1;
% ax.Position(2) = pos(2) - 0.027;
% ax.Position(3) = 0.1;
% ax.Position(4) = 0;
% plot(ax, [0,1], [0,0], 'LineStyle', '-', ...
%     'LineWidth', 1, ...
%     'Color', [0, 1, 0]);
% hold on;
% scatter(ax, 0.5, 0, 'Marker', '^', ...
%     'MarkerFaceColor', [0, 1, 0], ...
%     'MarkerEdgeColor', [0, 1, 0]);
% hold off;
% ax.XTickLabel = [];
% ax.YTickLabel = [];
% ax.XLim = [0, 1];
% ax.YLim = [-0.1, 0.1];

% txt = annotation(hf, 'textbox');
% txt.String = 'Good Trajectory';
% txt.Position(1) = ax.Position(1) + ax.Position(3) + 0.02;
% txt.Position(2) = ax.Position(2) - 0.022;
% txt.Position(3) = 0.3;
% txt.Position(4) = 0.05;
% txt.Interpreter = 'latex';
% txt.LineStyle = 'none';
% txt.FontSize = FONT_SIZE;

% pos = ax.Position;
% ax = axes(hf);
% ax.Position(1) = 0.1;
% ax.Position(2) = pos(2) - 0.05;
% ax.Position(3) = 0.1;
% ax.Position(4) = 0;
% plot(ax, [0,1], [0,0], 'LineStyle', '-', ...
%     'LineWidth', 1, ...
%     'Color', [1, 0, 0]);
% hold on;
% scatter(ax, 0.5, 0, 'Marker', 'x', ...
%     'MarkerFaceColor', [1, 0, 0], ...
%     'MarkerEdgeColor', [1, 0, 0], ...
%     'LineWidth', 2);
% hold off;
% ax.XTickLabel = [];
% ax.YTickLabel = [];
% ax.XLim = [0, 1];
% ax.YLim = [-0.1, 0.1];

% txt = annotation(hf, 'textbox');
% txt.String = 'Bad Trajectory';
% txt.Position(1) = ax.Position(1) + ax.Position(3) + 0.02;
% txt.Position(2) = ax.Position(2) - 0.022;
% txt.Position(3) = 0.3;
% txt.Position(4) = 0.05;
% txt.Interpreter = 'latex';
% txt.LineStyle = 'none';
% txt.FontSize = FONT_SIZE;

print(hf, '-r300', '-dpng', 'exampleFigs/pngs/cwh-example.png');

