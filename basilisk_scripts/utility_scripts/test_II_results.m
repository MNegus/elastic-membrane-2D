%% 
clear;
close all;

parentDir = "/media/michael/newarre/elastic_membrane/numbered_tests"
testNums = 1 : 2;

dirNames = strings(2, 1);
displayNames = strings(4, 1);
for num = testNums
    if (num == 1)
        dirNames(num) = sprintf("%s/centreUnadjustedAMR", parentDir);
        displayNames(num) = "Unadjusted";
    else
        dirNames(num) = sprintf("%s/test_II%d", parentDir, num);
        displayNames(num) = sprintf("Test II%d", num - 1);
    end
    
end

%% Exact solution
MEMBRANE_RADIUS = 3.0;
mag = 1.0;
WExact = @(x) mag * (1 - x.^2 / MEMBRANE_RADIUS^2);
xCentre = 0.25 * 6;
yCentre = 0.25 * 6;
% origCurve = @(x, y) (x - xCentre).^2 + (y - yCentre - WExact(x)).^2 - 1;
origCurve = @(x, y) (x - xCentre).^2 + (y - yCentre).^2 - 1;


%% Loop over all the test cases
for dirIdx = 1 : length(dirNames)
    outputFilename = sprintf("%s/raw_data/interface_10.txt", dirNames(dirIdx));
    
    % Load in the unsorted matrix A
    A = readmatrix(outputFilename);

    %% Scatter plot the interface shapes
    figure(1);
    hold on;
    scatter(A(:, 1), A(:, 2) - WExact(A(:, 1)), [], 'Displayname', displayNames(dirIdx));
    if (dirIdx == length(dirNames))
        fimplicit(origCurve, [0 5 0 5], 'linewidth', 2, 'color', 'black', ...
            'Linestyle', ':', 'Displayname', 't = 0');
    end

    
    ylim([1 4]);
    xlim([0 3]);
    pbaspect([1 1 1]);
    grid on;
    xlabel("$X$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$Y$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
    set(gcf, 'position', [200 200 800 800]);
    legend("interpreter", "latex", "Fontsize", 15);
    title("Tests II: Droplet in middle", "interpreter", "latex", "Fontsize", 18);

    %% Sort the matrix and determine the curvature



end