function setPublicationDefaults()
% setPublicationDefaults sets default plot properties for publication quality.
%
%   Call this function at the start of your session to apply the settings.
%
%   Example:
%       >> setPublicationDefaults();
%
%   This function resets previous defaults and sets:
%       - Figure background to white.
%       - Axes fonts to Times New Roman with size 14.
%       - Axes line width to 1.5 and tick direction to out.
%       - Line width to 2 and marker size to 8.
%       - A predefined color order.
%
%   Author: Your Name
%   Date:   Today's Date

% Reset to factory defaults
reset(groot);

% Set default figure properties.
set(groot, 'DefaultFigureColor', 'w');  % White background

% Set default axes properties.
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 30);
set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesBox', 'on');

% Set default line properties.
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 8);

% Set default text properties.
set(groot, 'DefaultTextFontName', 'Times New Roman');
set(groot, 'DefaultTextFontSize', 30);

% Set default legend properties.
set(groot, 'DefaultLegendFontSize', 26);

% Set a default color order for distinct plot colors.
set(groot, 'DefaultAxesColorOrder', [0, 0.4470, 0.7410;
                                     0.8500, 0.3250, 0.0980;
                                     0.9290, 0.6940, 0.1250;
                                     0.4940, 0.1840, 0.5560;
                                     0.4660, 0.6740, 0.1880;
                                     0.3010, 0.7450, 0.9330;
                                     0.6350, 0.0780, 0.1840]);

disp('Publication quality default settings have been applied.');
end
