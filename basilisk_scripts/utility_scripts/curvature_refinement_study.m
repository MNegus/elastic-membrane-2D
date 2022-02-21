%% curvature_refinement_study.m
% Validation of the curvature by altering the refinement level

clear;
close all;

% Radius and centre of the droplet in the lab frame
R = 1;
Yc = 0.125 + R;

% Grid definitions
boxWidth = 3;

% Membrane parameters
L = 1.25; % Width of membrane

% Range of levels
levels = 5 : 13;
 
% Set up figure
close all;
figure(1);
hold on;

%% Loop over levels
mkrs = ['o', '^', 's'];
mkrCtr = 1;
for shape = ["flat", "linear", "quadratic"]
    
    if shape == "flat"
        mag = 0;
    else
        mag = 0.5;
    end
    
    % Set W solution
%     WExact = @(X) mag * cos(pi * X / (2 * L)); 
%     WxExact = @(X) -(pi / (2 * L)) * mag * sin(pi * X / (2 * L));
%     WxxExact = @(X) -(pi / (2 * L))^2 * mag * cos(pi * X / (2 * L));
    WExact = @(X) mag * X.^2; 
    WxExact = @(X) 2 * mag * X;
    WxxExact = @(X) 2 * mag * ones(size(X));

    % L2 norm error arrays
    lineCentredL2Norms = zeros(size(levels));
    cellCentredL2Norms = zeros(size(levels));
    for levelIdx = 1 : length(levels)
%     for levelIdx = length(levels)
        level = levels(levelIdx);
        
        % Saves minimum cell size
        minCellSize = boxWidth / 2^level;

        % Load in output file
        outputFilename = sprintf("refinement_data_%s/interface_%d.txt", shape, level);
        A = readmatrix(outputFilename);

        % Set constants related to data size
        noSegments = size(A, 1) % Number of line segments
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
        
        % Exclude the curvature at the boundary
%         zeroIdxs = (xStarts == 0) | (xEnds == 0);
        zeroIdxs = (xStarts <= 3 * minCellSize) | (xEnds <= 3 * minCellSize);
        find(zeroIdxs)
        kappa(zeroIdxs) = nan;
        kappax(zeroIdxs) = nan;
        kappay(zeroIdxs) = nan;

%         close all;

        %% Plot curvature
%         subplot(2,1,1);
%         scatter(arcLengths, kappa, [], [0.9290, 0.6940, 0.1250]);
%         grid on
%         legend("kappa", "interpreter", "latex", "Fontsize", 15);
%         xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
%         ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
% %         ylim([2 - 1e-3, 2 + 1e-3]);
%         title("Comparison of curvature", "interpreter", "latex", "Fontsize", 18);
% 
%         subplot(2,1,2); 
%         scatter(arcLengths, kappax);
%         hold on;
%         scatter(arcLengths, kappay);
%         hold off;
%         grid on
%         legend(["kappa\_x", "kappa\_y"], "interpreter", "latex", "Fontsize", 15);
%         xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
%         ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
% %         ylim([2 - 1e-3, 2 + 1e-3]);
%         set(gcf, 'position', [200 200 1200 800]);
%         pause(2);

        %% Find centres of segments
        % Different to xCentres and yCentres, which are the x and y coordinates of
        % the centre of the CELL, not the segment
        xLineCentres = 0.5 * (xStarts + xEnds);
        yLineCentres = 0.5 * (yStarts + yEnds);

        %% Determine kappax in the cell centred or line centred formulation
        % Cell centred
        kappax_cellCentred = 2 * abs((hyy + hy.^3 .* WxxExact(xCellCentres)) ...
            ./ ((1 - WxExact(xCellCentres) .* hy).^2 + hy.^2).^(3/2));
        kappax_cellCentred(kappa == kappay) = nan;

        kappax_lineCentred = 2 * abs((hyy + hy.^3 .* WxxExact(xLineCentres)) ...
            ./ ((1 - WxExact(xLineCentres) .* hy).^2 + hy.^2).^(3/2));
        kappax_lineCentred(kappa == kappay) = nan;

        %% Determine L2-norm errors
        format long
%         cellCentredL2Norms(levelIdx) = sqrt(sum((kappax_cellCentred - 2).^2, 'omitnan') / noSegments);
        cellCentredL2Norms(levelIdx) = sqrt(sum((kappa - 2).^2, 'omitnan') / noSegments);
%         lineCentredL2Norms(levelIdx) = sqrt(sum((kappax_lineCentred - 2).^2, 'omitnan') / noSegments)

    end

    %% Plot L2-norm errors
    sz = 100;
    scatter(2.^levels, cellCentredL2Norms, sz, mkrs(mkrCtr), 'filled');
    mkrCtr = mkrCtr + 1;
%     plot(levels, lineCentredL2Norms, 'Displayname', sprintf("Line centred (%s)", shape));

end

%%
coeff2 = (2^levels(1))^2 * 0.044;
coeff1 = (2^levels(1)) * 0.168;
plot(2.^levels, coeff2 ./ (2.^levels).^2, 'linestyle', '--', ...
            'color', 'black', 'linewidth', 2, ...
            'Displayname', "$\sim 1 / N^2$");
plot(2.^levels, coeff1 ./ (2.^levels), 'linestyle', ':', ...
            'color', 'black', 'linewidth', 2, ...
            'Displayname', "$\sim 1 / N$");
        
% plot(levels, 28 * exp(-sqrt(2) * levels), 'linestyle', '--', ...
%             'color', 'black', 'linewidth', 2, ...
%             'Displayname', "$\sim \exp(-\sqrt{2} m)$)");
% plot(levels, 1.75 * exp(-levels / sqrt(2)), 'linestyle', ':', ...
%             'color', 'black', 'linewidth', 2, ...
%             'Displayname', "$\sim \exp(-m / \sqrt{2})$)");
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
grid on
legend(["Flat membrane", "Linear membrane", "Quadratic membrane", "$\sim 1 / N^2$", "$\sim 1 / N$"], "interpreter", "latex", "Fontsize", 15, 'location', 'southwest');
xlabel("Points per unit radius, $N$", "interpreter", "latex", "Fontsize", 18);
ylabel("RMS error", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
set(gcf, 'position', [200 200 800 600]);

exportgraphics(gcf, "curvature_figures/curvature_refinement.png");
savefig(gcf, "curvature_figures/curvature_refinement.fig");

