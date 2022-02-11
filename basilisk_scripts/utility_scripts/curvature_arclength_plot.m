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

%% Plot interface
figure(1);
plot(B(:, 1), B(:, 2));

%% Plot curvature as a function of arclength
kappa = B(:, 5);
kappax = B(:, 6);
kappay = B(:, 7);

% Determine indices where kappa is equal to kappax and kappay
xIdxs = (kappa == kappax);
yIdxs = (kappa == kappay);
% kappax(kappax > 1e3) = nan;
% kappay(kappay > 1e3) = nan;

close(figure(2));
figure(2);
hold on;
scatter(arcLengths(xIdxs), kappax(xIdxs));
scatter(arcLengths(yIdxs), kappay(yIdxs));
% plot(arcLengths, kappa, 'linewidth', 2);
legend(["kappax", "kappay", "kappa"]);


% Find non-NaN kappay values
nonNaNIdxs = not(isnan(kappay));

% Find non-2 values of kappay
nonTwoIdxs = not(kappay == 2);
idxs = nonNaNIdxs & nonTwoIdxs

% close(figure(3));
% figure(3);
% scatter(arcLengths(idxs), kappay(idxs), [], 'red');
set(gca, 'yscale', 'log');

%% Plot heights
hx = B(:, 8);
hxx = B(:, 9);
hy = B(:, 10);
hyy = B(:, 11);

% hx(abs(hx) > 1e3) = nan;
% hy(abs(hy) > 1e3) = nan;
% hxx(abs(hxx) > 1e3) = nan;
% hyy(abs(hyy) > 1e3) = nan;

close(figure(3));
figure(3);
hold on;
scatter(arcLengths(xIdxs), hy(xIdxs));
scatter(arcLengths(yIdxs), hx(yIdxs));
legend(["hy", "hx"]);

close(figure(4));
figure(4);
hold on;
scatter(arcLengths(xIdxs), hyy(xIdxs));
legend(["hyy", "hxx"]);
scatter(arcLengths(yIdxs), hxx(yIdxs));

%%
figure(5);
hold on;
scatter(arcLengths, xIdxs)
scatter(arcLengths, yIdxs)