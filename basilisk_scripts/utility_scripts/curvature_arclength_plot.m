%% curvature_arclength_plot.m
% Plots the curvature (amongst other quantities) of an interface as a
% function of arclength, so that multivalued interfaces can be represented.
% EXPLAIN THE INPUTS, METHOD AND OUTPUTS

clear;
close all;

% Load in output file
outputFilename = "interface_mag_0.75.txt";
A = readmatrix(outputFilename);

% Set constants related to data size
noSegments = size(A, 1); % Number of line segments
noVariables = size(A, 2) - 4; % Number of variables (e.g. curvature) 

% Initialise the matrix B, which will contain the data from A but sorted so
% that the line segments connect to each other
B = zeros(size(A));

%% Finding the first element
% Here we move the elements of A into B, sorting such that the line
% segments in consecutive rows are connected

% Parameters used in the algorithm
tol = 1e-2; % Tolerance to say two points are close (depends on level)
hugeVal = 1e5; % Huge value for distance calculations

% Find the index of the segment with the largest value of y (start or end
% point)
[maxStart, startIdx] = max(A(:, 2));
[maxEnd, endIdx] = max(A(:, 4));
if (maxStart > maxEnd)
    initIdx = startIdx;
else
    initIdx = endIdx;
end

% Copy the row from A over to B, and replace the row in A
B(1, :) = A(initIdx, :);
A(initIdx, :) = hugeVal * ones(1, noVariables + 4);

%% Sort the remaining rows 
startDists = zeros(noSegments); % Intialise distance array
currPoint = B(1, 3 : 4); % Current point 

for k = 2 : noSegments
% for k = 2
    k
   
    % End point of current segment
%     currPoint = B(k - 1, 3 : 4); 
    
    % Find distances of the start and end points in A form currPoint
    startDists = sqrt((currPoint(1) - A(:, 1)).^2 + (currPoint(2) - A(:, 2)).^2);
    endDists = sqrt((currPoint(1) - A(:, 3)).^2 + (currPoint(2) - A(:, 4)).^2);
    
    % Find the minimum distance
    [minStart, startIdx] = min(startDists);
    [minEnd, endIdx] = min(startDists);
    minVal = min([minStart, minEnd]);
    
    if (minStart < minEnd)
        % If the closest point is a start point, then update the new
        % currPoint to be the end point of the current segment
        idx = startIdx;
        currPoint = A(idx, 3 : 4);
    else
        % Else the closest point is an end point, so update the new
        % currPoint to be the start point of the current segment
        idx = endIdx;
        currPoint = A(idx, 1 : 2);
    end
%     
    
    % If minVal is larger than tol, then quit
    if (minVal > tol)
        disp("minVal = " + minVal);
        error("minVal is too large");
    end
    
    % Copy the row over to B and replace the row in A
    B(k, :) = A(idx, :);
    A(idx, :) = hugeVal * ones(1, noVariables + 4);
end


%% Determine arc length
segmentLengths = sqrt((B(:, 1) - B(:, 3)).^2 + (B(:, 2) - B(:, 4)).^2);
arcLengths = cumtrapz(segmentLengths);

%% Load in variables
kappa = B(:, 5);
kappax = B(:, 6);
kappay = B(:, 7);
heightY = B(:, 8);
hx = B(:, 9);
hxx = B(:, 10);
heightX = B(:, 11);
hy = B(:, 12);
hyy = B(:, 13);
W = B(:, 14);
Wx = B(:, 15);
Wxx = B(:, 16);

%% Plot interface
figure(1);
plot(B(:, 1), B(:, 2), 'linewidth', 2);
xlim([0 3]);
grid on;
xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
ylabel("$y$", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
title("Droplet interface in the moving frame", "interpreter", "latex", "Fontsize", 15);
pbaspect([1 1 1]); 
exportgraphics(gca, "curvature_figures/interface.png");
savefig(gcf, "curvature_figures/interface.fig");

%% Plot curvature as a function of arclength
close all;
% Determine indices where kappa is equal to kappax and kappay
kappax(kappax > 1e3) = nan;
kappay(kappay > 1e3) = nan;

subplot(2,1,1);
scatter(arcLengths, kappa, [], [0.9290, 0.6940, 0.1250]);
grid on
legend("kappa", "interpreter", "latex", "Fontsize", 15);
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
ylim([1.997, 2.003]);
title("Comparison of curvature", "interpreter", "latex", "Fontsize", 18);

subplot(2,1,2); 
hold on;
scatter(arcLengths, kappax);
scatter(arcLengths, kappay);
grid on
legend(["kappa\_x", "kappa\_y"], "interpreter", "latex", "Fontsize", 15);
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
ylim([1.997, 2.003]);


set(gcf, 'position', [200 200 1200 800]);

exportgraphics(gcf, "curvature_figures/kappas.png");
savefig(gcf, "curvature_figures/kappas.fig");


%% Plot height functions
heightY(abs(heightY) > 1e3) = nan;
heightX(abs(heightX) > 1e3) = nan;

close(figure(2));
figure(2)
hold on;
scatter(arcLengths, heightX);
scatter(arcLengths, heightY);
grid on
legend(["h.x[]", "h.y[]"], "interpreter", "latex", "Fontsize", 15);
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Height function", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);

title("Comparison of height functions", "interpreter", "latex", "Fontsize", 18);

set(gcf, 'position', [200 200 1400 600]);

exportgraphics(gcf, "curvature_figures/heights.png");
savefig(gcf, "curvature_figures/heights.fig");

%% Height function first derivatives
hx(abs(hx) > 1e3) = nan;
hy(abs(hy) > 1e3) = nan;

close(figure(3));
figure(3);
hold on;
scatter(arcLengths, hy);
scatter(arcLengths, hx);
grid on
legend(["hy (y derivative of h.x[])", "hx (x derivative of h.y[])"], ...
    "interpreter", "latex", "Fontsize", 15, 'location', 'northwest');
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Height derivative", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
set(gcf, 'position', [200 200 1400 400]);

title("First derivatives of height functions", "interpreter", "latex", "Fontsize", 18);


exportgraphics(gcf, "curvature_figures/height_first_derivative.png");
savefig(gcf, "curvature_figures/height_first_derivative.fig");


%% Height function second derivatives
hxx(abs(hxx) > 1e3) = nan;
hyy(abs(hyy) > 1e3) = nan;

close(figure(4));
figure(4);
hold on;
scatter(arcLengths, hyy);
scatter(arcLengths, hxx);
grid on
legend(["hyy (Second y derivative of h.x[])", "hxx (Second x derivative of h.y[])"], ...
    "interpreter", "latex", "Fontsize", 15, 'location', 'northwest');
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Height derivative", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
set(gcf, 'position', [200 200 1400 400]);

title("Second derivatives of height functions", "interpreter", "latex", "Fontsize", 18);

exportgraphics(gcf, "curvature_figures/height_second_derivative.png");
savefig(gcf, "curvature_figures/height_second_derivative.fig");


%% First derivative of W
close(figure(5));
figure(5);
scatter(arcLengths, Wx);
grid on
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Wx", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
set(gcf, 'position', [200 200 1400 400]);

title("First derivative of membrane position", "interpreter", "latex", "Fontsize", 18);

exportgraphics(gcf, "curvature_figures/membrane_first_derivative.png");
savefig(gcf, "curvature_figures/memrbane_first_derivative.fig");

%% Second derivative of W
close(figure(6));
figure(6);
scatter(arcLengths, Wxx);
grid on
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Wxx", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
set(gcf, 'position', [200 200 1400 400]);

title("Second derivative of membrane position", "interpreter", "latex", "Fontsize", 18);

exportgraphics(gcf, "curvature_figures/membrane_second_derivative.png");
savefig(gcf, "curvature_figures/memrbane_second_derivative.fig");

% %%
% % Compare kappax to its components
% 
% 
% kappax_manual = 2 * abs((hyy + hy.^3 .* Wxx) ./ ((1 - Wx .* hy).^2 + hy.^2).^(3/2));
% % kappax_manual = 2 * abs((hyy + hy.^3 .* Wxx) ./ ((1 - Wx).^2 + hy.^2).^(3/2));
% kappax_manual((abs(hy) > 1e3) | (abs(hyy) > 1e3)) = nan;
% 
% close(figure(6));
% figure(6)
% hold on;
% scatter(arcLengths, kappax);
% % scatter(arcLengths, kappax_manual);
% % scatter(arcLengths, Wx);
% % scatter(arcLengths, Wxx);
% % scatter(arcLengths, hy);
% scatter(arcLengths, hyy);
% % scatter(arcLengths, kappax_manual);
% % ylim([1.995, 2.005]);
% legend(["kappax", "kappax manual", "Wx", "Wxx", "hy", "hyy"]);
% 
% %%
% % Compare kappay to its components
% close(figure(7));
% figure(7)
% hold on;
% scatter(arcLengths, kappay);
% % scatter(arcLengths, Wx);
% % scatter(arcLengths, Wxx);
% scatter(arcLengths, hx);
% % scatter(arcLengths, hxx);
% % scatter(arcLengths, kappax_manual);
% legend(["kappay", "Wx", "Wxx", "hx", "hxx"]);

