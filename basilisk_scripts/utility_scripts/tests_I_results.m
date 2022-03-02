%% 
clear;
close all;

parentDir = "/media/michael/newarre/elastic_membrane/numbered_tests"
testNums = 1 : 5;

dirNames = strings(5, 1);
displayNames = strings(5, 1);
for num = testNums
    if (num == 1)
        dirNames(num) = sprintf("%s/wallUnadjustedAMR", parentDir);
        displayNames(num) = "Unadjusted";
    else
        dirNames(num) = sprintf("%s/test_I%d", parentDir, num - 1);
        displayNames(num) = sprintf("Test I%d", num - 1);
    end
    
end

%% Exact solution
MEMBRANE_RADIUS = 3.0;
mag = 1.0;
WExact = @(x) mag * (1 - x.^2 / MEMBRANE_RADIUS^2);
yCentre = 0.25 * 6;
origCurve = @(x, y) x.^2 + (y - yCentre - WExact(x)).^2 - 1;




%% Loop over all the test cases
for dirIdx = 1 : length(dirNames)
    outputFilename = sprintf("%s/raw_data/interface_10.txt", dirNames(dirIdx));
    
    % Load in the unsorted matrix A
    A = readmatrix(outputFilename);

    %% Scatter plot the interface shapes
    figure(1);
    hold on;
    scatter(A(:, 1), A(:, 2), [], 'Displayname', displayNames(dirIdx));
    if (dirIdx == 4)
        fimplicit(origCurve, [0 3 0 4], 'linewidth', 2, 'color', 'black', ...
            'Linestyle', ':', 'Displayname', 't = 0');
    end

    
    ylim([1 4]);
    xlim([0 3]);
    pbaspect([1 1 1]);
    grid on;
    xlabel("$x$", "interpreter", "latex", "Fontsize", 18);
    ylabel("$y$", "interpreter", "latex", "Fontsize", 18);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", 15);
    legend("interpreter", "latex", "Fontsize", 15);
    title("Tests I: Droplet at wall", "interpreter", "latex", "Fontsize", 18);
    set(gcf, 'position', [200 200 800 800]);

    %% Sort the matrix and determine the curvature



end