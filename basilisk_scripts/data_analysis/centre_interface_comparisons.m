%% interface_comparisons
close all;
clear;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex')

% Flag options
MOVING = 1;
MEMBRANE = 1;
WALL = 0;

position = "centre";



% parent_dir = sprintf("/media/michael/newarre/elastic_membrane/%s_tests_8_3_22", position);
parent_dir = "/home/negus/Desktop/cosTests";


% Data directories
untransformedDir = sprintf("%s/lab/code/spuriousTransformed", parent_dir);
transformedDir = sprintf("%s/curvilinear/code/spuriousTransformed", parent_dir);

% Variables
BOX_WIDTH = 6;
if WALL
    xCentre = 0;
    yCentre = 0.5 * BOX_WIDTH;
    L = 1.5;
    mag = 0.5;
else
    xCentre = 3.0;
    yCentre = 3.0;
    L = 5.0;
    mag = 0.5;
end


% Exact solution
% WExact = @(x) mag * (1 - x.^2 / L^2);
WExact = @(x) 0.5 * mag * (cos(pi * x / BOX_WIDTH) - 1);


for TIMESTEP = 0 : 191
%% Loads in untransformed case
outputFilename = sprintf("%s/interface_%d_9.txt", untransformedDir, TIMESTEP);
    
% Load in the unsorted matrix A
untransformedA = readmatrix(outputFilename);
XsUntransformed = untransformedA(:, 1);
YsUntransformed = untransformedA(:, 2);

%% Loads in transformed case
outputFilename = sprintf("%s/interface_%d_9.txt", transformedDir, TIMESTEP);
    
% Load in the unsorted matrix A
transformedA = readmatrix(outputFilename);
XsTransformed = transformedA(:, 1);
YsTransformed = transformedA(:, 2);

%% Plots lab frame
subplot(1, 2, 1);


% Plots untransformed curve 
scatter(XsUntransformed, YsUntransformed - WExact(XsUntransformed)); 
hold on;

% Plots transformed curve
scatter(XsTransformed, YsTransformed - WExact(XsTransformed)); 

% Plots exact solution
labCurve = @(x, y) (x - xCentre).^2 + (y - yCentre).^2 - 1;
fimplicit(labCurve, [0 BOX_WIDTH 0 BOX_WIDTH], 'linewidth', 2, 'color', 'black', ...
            'Linestyle', ':', 'Displayname', 'Target solution');
        
% Plot membrane
xs = linspace(0, L, 1e3);
WVals = WExact(xs);
scatter(xs, -WVals);
xline(0, 'color', 'black', 'linewidth', 2);
yline(0, 'color', 'black', 'linewidth', 2);

xlim([-max(WVals), BOX_WIDTH]);
ylim([-max(WVals), BOX_WIDTH]);

hold off;
        
grid on;
pbaspect([1 1 1]);
xticks(1:6);
xlabel('$x$');
ylabel('$y$');

legend('Untransformed equations', 'Transformed equations', 'Target solution');
title("Lab frame");

%% Plots curvilinear frame
subplot(1, 2, 2);

% Plots untransformed curve 
scatter(XsUntransformed, YsUntransformed); 
hold on;

% Plots transformed curve
scatter(XsTransformed, YsTransformed); 

% Plots exact solution
curvilinearCurve = @(x, y) (x - xCentre).^2 + (y - yCentre - WExact(x)).^2 - 1;
fimplicit(curvilinearCurve, [0 BOX_WIDTH 0 BOX_WIDTH], 'linewidth', 2, 'color', 'black', ...
            'Linestyle', ':', 'Displayname', 'Target solution');
        
xlim([-max(WVals), BOX_WIDTH]);
ylim([-max(WVals), BOX_WIDTH]);
scatter(xs, zeros(size(xs)));
xline(0, 'color', 'black', 'linewidth', 2);
yline(0, 'color', 'black', 'linewidth', 2);

hold off;

grid on;
pbaspect([1 1 1]);
xticks(1:6);
xlabel('$X$');
ylabel('$Y$');

legend('Untransformed equations', 'Transformed equations', 'Target solution');
title("Curvilinear frame");

%% General figure properties
sgtitle("$t$ =" + num2str(TIMESTEP) , 'Fontsize', 20);
set(gcf, 'position', [200 200 1200 600]);

drawnow;
pause(0.01);
TIMESTEP

end
