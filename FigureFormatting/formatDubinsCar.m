% 
% Name: formatDubinsCar.m
% Author: Joseph D. Gleason
% 
% Description: Format the Dubins car example figure for paper plotting
% 

close all;
clearvars;

% paper column width in pts
HSCC_COLWIDTH = 241.14749;

FONT_SIZE = 8;

% Load the figure
hf = openfig('../exampleFigs/DubinsCar_example_point.fig');
delete(hf.Children(1));
ha = hf.Children(1);

% undock the figures
hf.WindowStyle = 'normal';

% need to pause to let MATLAB catch up to the undocking
pause(2);

% reposition figures
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 2.8;
hf.Position(1) = 3;
hf.Position(2) = 3;

% get rid of ratioed scaling
ha.DataAspectRatioMode = 'auto';
ha.PlotBoxAspectRatioMode = 'auto';

% Reset axes details
ha.FontSize = FONT_SIZE;
ha.TickLabelInterpreter = 'latex';
ha.XLabel.String = '$x$';
ha.XLabel.FontSize = FONT_SIZE;
ha.XLabel.Interpreter = 'latex';

ha.YLabel.String = '$y$';
ha.YLabel.FontSize = FONT_SIZE;
ha.YLabel.Interpreter = 'latex';

delete(ha.Children(6));

scat_colors = {[1, 0, 0], [0, 0, 1], [0, 0.9, 0], [0, 0, 0], [1, 0.6, 0]};

scats = 0;
for lv = 1:length(ha.Children)
    if strcmp(ha.Children(lv).Type, 'patch')
        ha.Children(lv).LineStyle = '-';
        ha.Children(lv).EdgeColor = 0.7 * ones(1, 3);
        ha.Children(lv).EdgeAlpha = 0.5;
    end

    if strcmp(ha.Children(lv).Type, 'scatter')
        scats = scats + 1;
        ha.Children(lv).MarkerEdgeColor = 'none';
        ha.Children(lv).MarkerFaceColor = scat_colors{scats};
    end
end

ha.Children(length(ha.Children)).LineStyle = '-';
ha.Children(length(ha.Children)).EdgeColor = zeros(1, 3);
ha.Children(length(ha.Children)).EdgeAlpha = 1;

% ha.XLim = [-7.5, 9.3];
ha.Units = 'inches';
ha.Position(4) = 2.5;
ha.Position(3) = 1.8;
ha.Position(3:4) = ha.Position(3:4) * 0.9;
ha.Units = 'normalized';
ha.Position(2) = 0.13;
ha.Position(1) = 0.11;
% ha.Position(2) = 0.45;

hl = legend([ha.Children(1:5)], ...
    {'\texttt{lag-under}' '\texttt{particle-open}', '\texttt{genzps-open}', ...
        '\texttt{chance-affine}', '\texttt{chance-open}'});

hl.FontSize = FONT_SIZE;
hl.Interpreter = 'latex';
hl.Box = 'off';
hl.Position = [0.5803    0.5984    0.3845    0.2108];

an = annotation(hf, 'rectangle');
an.Position(1) = hl.Position(1) + 0.062;
an.Position(2) = hl.Position(2) - 0.07;
an.Units = 'inches';
an.Position(3) = 0.1;
an.Position(4) = 0.1;
an.Units = 'normalized';
an.FaceColor = 'y';
pos = an.Position;

for lv = 1:5
    antemp = annotation(hf, 'rectangle');
    antemp.Position(1) = pos(1) + 0.0008;
    antemp.Position(2) = pos(2) - 0.005;
    % antemp.Units = 'inches';
    antemp.Position(3) = pos(3)/1.08;
    antemp.Position(4) = pos(4)/1.08;
    % antemp.Units = 'normalized';
    antemp.FaceColor = 'y';

    antemp.LineStyle = '-';
    antemp.EdgeColor = 0.7 * ones(1, 3);
    % antemp.EdgeAlpha = 0.5;
    antemp.FaceAlpha = 0.08;

    pos = antemp.Position;
end

an = annotation(hf, 'rectangle');
an.Position(1) = hl.Position(1) + 0.062;
an.Position(2) = hl.Position(2) - 0.07;
an.Units = 'inches';
an.Position(3) = 0.1;
an.Position(4) = 0.1;
an.Units = 'normalized';
an.FaceColor = 'y';
pos = an.Position;

txt = annotation(hf, 'textbox');
txt.LineStyle = 'none';
txt.Position(1) = an.Position(1) + 0.064;
txt.Position(2) = an.Position(2) + 0.001;
txt.Position(3) = 0.3;
txt.Position(4) = 0.05;
txt.String = 'Target Tube';
txt.Interpreter = 'latex';
txt.FontSize = FONT_SIZE;

print(hf, '-r300', '-dpng', '../exampleFigs/pngs/dubinscar-example.png');