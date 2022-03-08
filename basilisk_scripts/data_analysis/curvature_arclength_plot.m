%% curvature_arclength_plot.m
% Plots the curvature (amongst other quantities) of an interface as a
% function of arclength, so that multivalued interfaces can be represented.
% EXPLAIN THE INPUTS, METHOD AND OUTPUTS

clear;
close all;

% Radius and centre of the droplet in the lab frame
R = 1;
Yc = 0.125 + R;

% Membrane parameters
mag = 0.5; % Magnitude of displacement
L = 1.25; % Width of membrane

% Exact solution for membrane
% WExact = @(X) mag * cos(pi * X / (2 * L)); 
% WxExact = @(X) -(pi / (2 * L)) * mag * sin(pi * X / (2 * L));
% WxxExact = @(X) -(pi / (2 * L))^2 * mag * cos(pi * X / (2 * L));
WExact = @(X) mag * X.^2; 
WxExact = @(X) 2 * mag * X;
WxxExact = @(X) 2 * mag * ones(size(X));


% Load in output file
outputFilename = "interface_0.txt";
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
%     if (minVal > tol)
%         disp("minVal = " + minVal);
%         error("minVal is too large");
%     end
    
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

%% Find centres of segments
% Different to xCentres and yCentres, which are the x and y coordinates of
% the centre of the CELL, not the segment
xLineCentres = 0.5 * (xStarts + xEnds);
yLineCentres = 0.5 * (yStarts + yEnds);


%% Find exact solution for the interface, Y = T(X)
% As Y = T(X) will be multivalued, we define TPlus and TMinus to handle the
% upper and lower part of the droplet, respectively 

% Use the xCentres array
% X = xCentres;
% 
% TPlus = Yc + WExact(X) + sqrt(R^2 - X.^2);
% TMinus = Yc + WExact(X) - sqrt(R^2 - X.^2);
% 
% TPlusX = WxExact(X) - X ./ sqrt(R^2 - X.^2);
% TMinusX = WxExact(X) + X ./ sqrt(R^2 - X.^2);
% 
% TPlusXX = -WxxExact(X) - (1 ./ sqrt(R^2 - X.^2) ...
%     + X.^2 ./ (R^2 - X.^2).^(3/2));
% TMinusXX = -WxxExact(X) + (1 ./ sqrt(R^2 - X.^2) ...
%     + X.^2 ./ (R^2 - X.^2).^(3/2));
% 


%% Find exact solution for the interface, X = S(Y)
% Define range for Y values
% Ymin = Yc + WExact(0) - R;
% Ymax = Yc + WExact(0) + R;
% Y = linspace(Ymin, Ymax, 1e3);
% 
% % Use minimimsation to find X and a function of Y
% zeroFun = @(X) X.^2 + (Y - WExact(X) - Yc).^2 - R^2; % Levelset function
% X0 = sqrt(R^2 - (Y - Yc).^2); % Initial guess for X values
% S = fsolve(zeroFun, X0); % Solve for S, i.e. X = S(Y)
% S = abs(S); % Ensure we keep the positive solutions
% 
% % Find derivatives of S
% dy = Y(2) - Y(1);
% 
% % First Y derivative
% SY = zeros(size(S));
% SY(2 : end - 1) = (S(3 : end) - S(1 : end - 2)) / (2 * dy);
% SY(1) = (S(2) - S(1)) / dy;
% SY(end) = (S(end) - S(end - 1)) / dy;
% 
% % Second y derivative
% SYY = zeros(size(S));
% SYY(2 : end - 1) = (S(3 : end) - 2 * S(2 : end - 1) + S(1 : end - 2)) / dy^2;
% SYY(1) = (SY(2) - SY(1)) / dy;
% SYY(end) = (SY(end) - SY(end - 1)) / dy;

%% Plot kappa and arc lengths
% close all;
% figure(56565);
% hold on;
% scatter(yCentres, kappay);
% scatter(yCentres, kappax);
% % scatter(yCentres, kappa);


%% Plot interface
close all;
figure(1);
hold on;
% plot(B(:, 1), B(:, 2), 'linewidth', 2);
scatter(B(:, 1), B(:, 2));

% Plot exact interface
% plot(S, Y, 'color', 'black');

xlim([0 3]);
grid on;
xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
ylabel("$y$", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
title("Droplet interface in the moving frame", "interpreter", "latex", "Fontsize", 15);
pbaspect([1 1 1]); 
exportgraphics(gca, "curvature_figures/interface.png");
savefig(gcf, "curvature_figures/interface.fig");

%% Errors in hx

% % First isolate hx and hxx from its positive and negative parts
% poshxIdxs = (hx > 0) & (kappa == kappay) & (kappax ~= kappay);
% neghxIdxs = (hx < 0) & (kappa == kappay) & (kappax ~= kappay);
% 
% % Restrict numerical solutions for positive solution
% xPosq = xCentres(poshxIdxs);
% xNegq = xCentres(neghxIdxs);
% hxPosq = hx(poshxIdxs);
% hxNegq = hx(neghxIdxs);
% 
% 
% % Restrict the exact solution along xq
% TPlusXq = TPlusX(neghxIdxs);
% TMinusXq = TMinusX(poshxIdxs);
% 
% % Plot error
% close(figure(22));
% figure(22);
% hold on;
% scatter([xPosq; xNegq], [TMinusXq - hxPosq; TPlusXq - hxNegq]);
% grid on
% xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
% ylabel("Error", "interpreter", "latex", "Fontsize", 18);
% set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
% title("Error between $h_x$ and exact", "interpreter", "latex", "Fontsize", 18);
% set(gcf, 'position', [200 200 1400 600]);

%% Errors in hxx

% % First isolat hxx from its positive and negative parts
% poshxxIdxs = (hxx > 0) & (kappa == kappay) & (kappax ~= kappay);
% neghxxIdxs = (hxx < 0) & (kappa == kappay) & (kappax ~= kappay);
% 
% % Restrict numerical solutions for positive solution
% xPosq = xCentres(poshxxIdxs);
% xNegq = xCentres(neghxxIdxs);
% hxxPosq = hxx(poshxxIdxs);
% hxxNegq = hxx(neghxxIdxs);
% 
% 
% % Restrict the exact solution along xq
% TPlusXXq = TPlusXX(neghxxIdxs);
% TMinusXXq = TMinusXX(poshxxIdxs);
% 
% % Plot error
% close(figure(23));
% figure(23);
% hold on;
% % scatter([xPosq; xNegq], [TMinusXXq - hxxPosq; TPlusXXq - hxxNegq]);
% scatter([xPosq; xNegq], [TMinusXXq; TPlusXXq]);
% scatter([xPosq; xNegq], [hxxPosq; hxxNegq]);
% grid on
% xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
% ylabel("Error", "interpreter", "latex", "Fontsize", 18);
% set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
% title("Error between $h_{xx}$ and exact", "interpreter", "latex", "Fontsize", 18);
% set(gcf, 'position', [200 200 1400 600]);
% 

%% Errors in hy and hyy

% % Find where hy and hyy are nonzero and non-NaN, and such that the
% % contribution to kappa is from kappax.
% nonZeroYIdxs = (hy ~= 0) & (hyy ~= 0) & not(isnan(hy)) & not(isnan(hyy)) ...
%     & (kappa == kappax) & (kappax ~= kappay);
% 
% % Restrict the exact solutions to where hy and hy ~= 0
% yq = yCentres(nonZeroYIdxs);
% hyq = hy(nonZeroYIdxs);
% hyyq = hyy(nonZeroYIdxs);
% 
% % Plot errors between the exact solutions for height functions
% Syq = interp1(Y, SY, yq);
% Syyq = interp1(Y, SYY, yq);
% 
% close(figure(15));Line
% figure(15);
% plot(yq, abs(Syq - hyq));
% 
% close(figure(16));
% figure(16);
% plot(yq, abs(Syyq - hyyq));


%% Plot curvature as a function of arclength
close all;

subplot(2,1,1);
scatter(arcLengths, kappa, [], [0.9290, 0.6940, 0.1250]);
grid on
legend("kappa", "interpreter", "latex", "Fontsize", 15);
xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
ylim([0.29215, 0.29235]);
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
ylim([0.29215, 0.29235]);


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
hxx(abs(hxx) > 1e3) = nan;
hyy(abs(hyy) > 1e3) = nan;
title("Comparison of height functions", "interpreter", "latex", "Fontsize", 18);

set(gcf, 'position', [200 200 1400 600]);

exportgraphics(gcf, "curvature_figures/heights.png");
savefig(gcf, "curvature_figures/heights.fig");

%% Height function first derivatives

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

%%
% Compare kappax to its components
% kappax_cellCentred = 2 * abs((hyy + hy.^3 .* WxxExact(xCellCentres)) ./ ((1 - WxExact(xCellCentres) .* hy).^2 + hy.^2).^(3/2));
% % kappax_cellCentred(kappax_cellCentred == 0) = nan;
% kappax_cellCentred(kappa == kappay) = nan;
% 
% kappax_lineCentred = 2 * abs((hyy + hy.^3 .* WxxExact(xLineCentres)) ./ ((1 - WxExact(xLineCentres) .* hy).^2 + hy.^2).^(3/2));
% % kappax_lineCentred(kappax_lineCentred == 0) = nan;
% kappax_lineCentred(kappa == kappay) = nan;
% 
% 
% 
% close(figure(6));
% figure(6)
% hold on;
% scatter(arcLengths, kappax_cellCentred);
% scatter(arcLengths, kappax_lineCentred);
% % scatter(arcLengths, kappax_manual);
% grid on
% legend(["kappa (Cell centred)", "kappa (Line centred)"], ...
%     "interpreter", "latex", "Fontsize", 15, 'location', 'northeast');
% xlabel("Arc length", "interpreter", "latex", "Fontsize", 18);
% ylabel("Curvature", "interpreter", "latex", "Fontsize", 18);
% set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
% set(gcf, 'position', [200 200 800 600]);
% 
% % title("Second derivative of membrane position", "interpreter", "latex", "Fontsize", 18);
% exportgraphics(gcf, "curvature_figures/cell_vs_line_centred.png");
% savefig(gcf, "curvature_figures/cell_vs_line_centred.fig");

%% 
% close all;
% figure(45);
% hold on;
% originalTerm = hyy;
% originalTerm(kappa == kappay) = nan;
% adjustedTerm = hy.^3 .* Wxx;
% adjustedTerm(kappa == kappay) = nan;
% 
% % scatter(arcLengths, originalTerm);
% % scatter(arcLengths, adjustedTerm);
% 
% numerator = hyy + hy.^3 .* Wxx;
% numerator(kappa == kappay) = nan;
% 
% denominator = ((1 - Wx .* hy).^2 + hy.^2).^(3/2);
% denominator(kappa == kappay) = nan;
% 
% fullsol = numerator ./ denominator;
% fullsol(kappa == kappay) = nan;
% 
% % scatter(arcLengths, numerator);
% % scatter(arcLengths, denominator);
% scatter(arcLengths, numerator + denominator);
% 
% %% L2 norm error for kappax
% % Determines the L2-norm error for the kappax calulcations, given that it
% % should be 2 everywhere
% kappax_cellCentred_L2Norm = sqrt(sum((kappax_cellCentred - 2).^2, 'omitnan'))
% kappax_lineCentred_L2Norm = sqrt(sum((kappax_lineCentred - 2).^2, 'omitnan'))
% 
% 
% %%
% % Compare kappay to its components
% % close(figure(7));
% % figure(7)
% % hold on;
% % scatter(arcLengths, kappay);
% % % scatter(arcLengths, Wx);
% % % scatter(arcLengths, Wxx);
% % scatter(arcLengths, hx);
% % % scatter(arcLengths, hxx);
% % % scatter(arcLengths, kappax_manual);
% % legend(["kappay", "Wx", "Wxx", "hx", "hxx"]);
% 
