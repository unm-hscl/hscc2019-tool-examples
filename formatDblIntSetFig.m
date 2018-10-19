%
% Name:        formatDblIntSetFig.m
% Author:      Joseph D. Gleason
% 
% Descritpion: Format the double integrator set figure for publication
% 

close all;
clearvars;

% Set computation times
SET_COMP_TIMES = struct('lagunder', 0.22614, ...
                        'lagover', 0.21797, ...
                        'dynprog', 126.15970, ...
                        'ccc', 13.62379, ...
                        'ftgenz', 338.5651);

% paper column width in pts
HSCC_COLWIDTH = 241.14749;

FONT_SIZE = 8;

% Load the figure
openfig('exampleFigs/dblint_sets.fig');

% handles -> variables
hf = gcf;
ha = gca;
patches = ha.Children;

% Main figure setup
hf.Units = 'points';
hf.Position(3) = HSCC_COLWIDTH;
hf.Units = 'inches';
hf.Position(4) = 3.4;

% Axes
ha.Position(3) = 0.6;
ha.Units = 'inches';
ha.Position(4) = ha.Position(3);
ha.Units = 'normalized';
ha.Position(1) = 0.23;
ha.Position(2) = 0.37;
ha.XLim = [-1, 1];
ha.YLim = [-1, 1];

ha.FontSize = FONT_SIZE;
ha.TickLabelInterpreter = 'latex';
ha.XLabel.FontSize = FONT_SIZE;
ha.YLabel.FontSize = FONT_SIZE;

% the "Legend"
% Safe set
% -----------------------
% the box
an = annotation(hf, 'rectangle');
an.Position(1) = 0.1;
an.Position(2) = ha.Position(2) - 0.15;
an.Units = 'inches';
an.Position(3) = 0.1;
an.Position(4) = 0.1;
an.Units = 'normalized';
an.FaceColor = patches(6).FaceColor;

% Text
txt = annotation(hf, 'textbox');
txt.Position(1) = an.Position(1) + 0.03;
txt.Position(2) = an.Position(2) - 0.005;
txt.Position(3) = 0.3;
txt.Position(4) = 0.05;
txt.String = 'Safe Set';
txt.Interpreter = 'latex';
txt.LineStyle = 'none';
txt.FontSize = FONT_SIZE;

% comptimes text
cptxt = annotation(hf, 'textbox');
cptxt.Position(1) = txt.Position(1) + txt.Position(3) + 0.05;
cptxt.Position(2) = txt.Position(2);
cptxt.Position(3) = 0.5;
cptxt.Position(4) = 0.05;
cptxt.String = 'Computation Times [s]';
cptxt.Interpreter = 'latex';
cptxt.LineStyle = 'none';
cptxt.FontSize = FONT_SIZE;

% for loop the rest
strs = {'lagover', 'dynprog', 'ccc', 'ftgenz', 'lagunder'};
for lv = 1:length(strs)
    % the box
    pos = an.Position;
    an = annotation(hf, 'rectangle');
    an.Position(1) = pos(1);
    an.Position(2) = pos(2) - 0.04;
    an.Position(3) = pos(3);
    an.Position(4) = pos(4);
    an.Units = 'normalized';
    if lv < length(strs)
        an.FaceColor = patches(6-lv).FaceColor;
    else
        an.FaceColor = patches(6-lv+1).FaceColor;
        an = annotation(hf, 'rectangle');
        an.Position(1) = pos(1);
        an.Position(2) = pos(2) - 0.04;
        an.Position(3) = pos(3);
        an.Position(4) = pos(4);
        an.Units = 'normalized';
        an.FaceColor = patches(6-lv).FaceColor;
        an.FaceAlpha = 0.5;
    end

    % Text
    txt = annotation(hf, 'textbox');
    txt.Position(1) = an.Position(1) + 0.03;
    txt.Position(2) = an.Position(2) - 0.005;
    txt.Position(3) = 0.3;
    txt.Position(4) = 0.05;
    switch(strs{lv})
        case 'lagover'
            txt.String = '\texttt{lag-over}';
        case 'lagunder'
            txt.String = '\texttt{lag-under}';
        case 'dynprog'
            txt.String = '\texttt{SReachDynProg}';
        case 'ccc'
            txt.String = '\texttt{chance-open}';
        case 'ftgenz'
            txt.String = '\texttt{genzps-open}';
    end 
    txt.Interpreter = 'latex';
    txt.LineStyle = 'none';
    txt.FontSize = FONT_SIZE;

    % computation time
    cptxt = annotation(hf, 'textbox');
    cptxt.Position(1) = txt.Position(1) + txt.Position(3) + 0.05;
    cptxt.Position(2) = txt.Position(2);
    cptxt.Position(3) = 0.5;
    cptxt.Position(4) = 0.05;
    cptxt.String = num2str(SET_COMP_TIMES.(strs{lv}));
    cptxt.Interpreter = 'latex';
    cptxt.LineStyle = 'none';
    cptxt.FontSize = FONT_SIZE;
end

% % the box
% pos = an.Position;
% an = annotation(hf, 'rectangle');
% an.Position(1) = pos(1);
% an.Position(2) = pos(2) - 0.05;
% an.Position(3) = pos(3);
% an.Position(4) = pos(4);
% an.Units = 'normalized';
% an.FaceColor = patches(5).FaceColor;
% 
% % Text
% txt = annotation(hf, 'textbox');
% txt.Position(1) = an.Position(1) + 0.05;
% txt.Position(2) = an.Position(2) - 0.013;
% txt.Position(3) = 0.3;
% txt.Position(4) = 0.05;
% txt.String = '\texttt{lag-over}';
% txt.Interpreter = 'latex';
% txt.LineStyle = 'none';
% txt.FontSize = FONT_SIZE;
% 
% % computation time
% cptxt = annotation(hf, 'textbox');
% cptxt.Position(1) = txt.Position(1) + txt.Position(3) + 0.05;
% cptxt.Position(2) = txt.Position(2);
% cptxt.Position(3) = 0.5;
% cptxt.Position(4) = 0.05;
% cptxt.String = num2str(SET_COMP_TIMES.lagover);
% cptxt.Interpreter = 'latex';
% cptxt.LineStyle = 'none';
% cptxt.FontSize = FONT_SIZE;
% 
% % Dynamic Programming
% % ----------------------------------------
% % the box
% pos = an.Position;
% an = annotation(hf, 'rectangle');
% an.Position(1) = pos(1);
% an.Position(2) = pos(2) - 0.05;
% an.Position(3) = pos(3);
% an.Position(4) = pos(4);
% an.Units = 'normalized';
% an.FaceColor = patches(4).FaceColor;
% 
% % Text
% txt = annotation(hf, 'textbox');
% txt.Position(1) = an.Position(1) + 0.05;
% txt.Position(2) = an.Position(2) - 0.013;
% txt.Position(3) = 0.3;
% txt.Position(4) = 0.05;
% txt.String = '\texttt{SReachDynProg}';
% txt.Interpreter = 'latex';
% txt.LineStyle = 'none';
% txt.FontSize = FONT_SIZE;
% 
% % computation time
% cptxt = annotation(hf, 'textbox');
% cptxt.Position(1) = txt.Position(1) + txt.Position(3) + 0.05;
% cptxt.Position(2) = txt.Position(2);
% cptxt.Position(3) = 0.5;
% cptxt.Position(4) = 0.05;
% cptxt.String = num2str(SET_COMP_TIMES.dynprog);
% cptxt.Interpreter = 'latex';
% cptxt.LineStyle = 'none';
% cptxt.FontSize = FONT_SIZE;


print(hf, '-r300', '-dpng', 'exampleFigs/pngs/dlbint_sets.png');


