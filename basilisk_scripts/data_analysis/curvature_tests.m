%% curvature_refinement_study.m
% Validation of the curvature by altering the refinement level

clear;
close all;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex')

%% Parameters
% Parent directory
parent_dir = "/media/michael/newarre/elastic_membrane/curvatureRefineTests";

% Flag options
MOVING = 1;
MEMBRANE = 1;
WALL = 0;
TRANSPOSED = 1;

% Set keywords depending on parameters
if MOVING
    frame = 'Curvilinear';
else
    frame = 'Lab'
end

if WALL
    position = 'wall';
else
    position = 'embedded';
end

if TRANSPOSED
    transp = 'transposed';
else
    transp = 'original';
end

% Data directory
data_dir = sprintf("%s/MOVING_%d-MEMBRANE_%d-WALL_%d-TRANSPOSED_%d/raw_data", ...
    parent_dir, MOVING, MEMBRANE, WALL, TRANSPOSED)

% Grid definitions
BOX_WIDTH = 6;

% Flag dependent variables
if WALL
    if TRANSPOSED
        xCentre = 0.5 * BOX_WIDTH;
        yCentre = 0;
    else
        xCentre = 0;
        yCentre = 0.5 * BOX_WIDTH;
    end
    L = 1.5;
    
    mag = 0.5;
else
    xCentre = 3.0;
    yCentre = 3.0;
    L = 5.0;
    
    mag = 0.5;
end

% Range of levels
levels = 5 : 13;


% Exact solution for W
WExact = @(X) mag * (1 - X.^2 / L^2); 
WxExact = @(X) mag * (-2 * X / L^2);
WxxExact = @(X) mag * (-2 / L^2) * ones(size(X));


%% Loop over levels

% L2 norm error arrays
kappaL2Norms = zeros(size(levels));
kappaxL2Norms = zeros(size(levels));
kappayL2Norms = zeros(size(levels));

for levelIdx = 1 : length(levels)
    level = levels(levelIdx);

    % Saves minimum cell size
    minCellSize = BOX_WIDTH / 2^level;

    % Load in output file
    outputFilename = sprintf("%s/interface_0_%d.txt", data_dir, level);
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
    tol = 1e6; % Tolerance to say two points are close (depends on level)
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
        k

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
    xStarts = B(:, 1);
    yStarts = B(:, 2);
    xEnds = B(:, 3);
    yEnds = B(:, 4);
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
    xCellCentres = B(:, 17);
    yCellCentres = B(:, 18);

    % Convert large values to NaN
    kappax(kappax > 1e3) = nan;
    kappay(kappay > 1e3) = nan;
    hx(abs(hx) > 1e3) = nan;
    hy(abs(hy) > 1e3) = nan;
    hxx(abs(hxx) > 1e3) = nan;
    hyy(abs(hyy) > 1e3) = nan;

    %% Plot interface
    if (level == max(levels))
        
        subplot(1, 2, 1);
        hold on;
        h(1) = scatter(B(:, 1), B(:, 2));
        
        if (TRANSPOSED)
             ys = linspace(0, L, 1e3);
             WVals = WExact(ys);
             if MOVING
                 WVals = 0 * WVals;
             end
             h(2) = scatter(-WVals, ys);
             
             xlim([-max(WVals), BOX_WIDTH]);
             ylim([-max(WVals), BOX_WIDTH]);
        else
            xs = linspace(0, L, 1e3);
            WVals = WExact(xs);
            if MOVING
                 WVals = 0 * WVals;
            end
            h(2) = scatter(xs, -WVals);
            
            xlim([-max(WVals), BOX_WIDTH]);
            ylim([-max(WVals), BOX_WIDTH]);
        end
        
        h(3) = xline(0, 'color', 'black', 'linewidth', 2);
        h(4) = yline(0, 'color', 'black', 'linewidth', 2);
       
        if (MOVING)         
            % Axes labels
            xlabel('$X$');
            ylabel('$Y$');
        else
            xlabel('$x$');
            ylabel('$y$');
        end
        
        grid on;
        pbaspect([1 1 1]);
        legend(h(1:2), 'Droplet interface', 'Membrane position');        
        xticks(1:6);
        
%         figure(3);
% %         scatter(arcLengths, kappax);
%         scatter(arcLengths(kappa == kappax), kappax(kappa == kappax));
%         hold on;
% %         scatter(arcLengths, kappay);
%         scatter(arcLengths(kappa == kappay), kappay(kappa == kappay));
%         hold off;
%         grid on
%         legend(["kappa\_x", "kappa\_y"], "interpreter", "latex", "Fontsize", 15);
%         xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
%         ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
%         set(gcf, 'position', [200 200 1200 800]);
    end

    %% Determine L2-norm errors
    format long
    
    % kappa L2-norm
    kappaL2Norms(levelIdx) = sqrt(sum((kappa - 1).^2, 'omitnan')) / noSegments;
    
    % kappax L2-norm, only compare where kappax == kappa
    kappaxIdxs = kappax == kappa;
    kappaxVals = kappax(kappaxIdxs);
    nokappaxSegments = sum(~isnan(kappaxIdxs));
    kappaxL2Norms(levelIdx) = sqrt(sum((kappaxVals - 1).^2, 'omitnan')) / nokappaxSegments;
    
    % kappay L2-norm, only compare where kappay == kappa
    kappayIdxs = kappay == kappa;
    kappayVals = kappay(kappayIdxs);
    nokappaySegments = sum(~isnan(kappayVals));
    kappayL2Norms(levelIdx) = sqrt(sum((kappayVals - 1).^2, 'omitnan')) / nokappaySegments;
    

end

%% Plot L2-norm errors
subplot(1, 2, 2);
hold on;
sz = 100;
scatter(2.^levels, kappaL2Norms, sz, 'Filled', 'Displayname', 'kappa');
scatter(2.^levels, kappaxL2Norms, sz, 'd', 'Filled', 'Displayname', 'kappa\_x');
scatter(2.^levels, kappayL2Norms, sz, '^', 'Filled', 'Displayname', 'kappa\_y');

%% Err plotting
sampleIdx = 4;

if TRANSPOSED
    coeff2 = (2^levels(sampleIdx))^2 * kappaxL2Norms(sampleIdx);
    coeff1 = (2^levels(sampleIdx)) * kappayL2Norms(sampleIdx);
else
    coeff2 = (2^levels(sampleIdx))^2 * kappayL2Norms(sampleIdx);
    coeff1 = (2^levels(sampleIdx)) * kappaxL2Norms(sampleIdx);
end


plot(2.^levels, coeff1 ./ (2.^levels), 'linestyle', '--', ...
            'color', 'black', 'linewidth', 2, ...
            'Displayname', "$\sim 1 / N$");
plot(2.^levels, coeff2 ./ (2.^levels).^2, 'linestyle', ':', ...
            'color', 'black', 'linewidth', 2, ...
            'Displayname', "$\sim 1 / N^2$");

legend('Location', 'Southwest');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
grid on;
xlabel("Points per unit radius, $N$");
ylabel("RMS error of curvature");


% Set title dependent on parameters
titlestr = sprintf("%s frame, %s droplet, %s coordinates", ...
    frame, position, transp);
sgtitle(titlestr, 'Fontsize', 20); 

set(gcf, 'position', [200 200 1200 600]);

figName = sprintf("curvatureTests/%s_frame_%s_droplet_%s_coordinates", ...
    frame, position, transp);

exportgraphics(gcf, sprintf("%s.png", figName));
savefig(gcf, sprintf("%s.fig", figName));

