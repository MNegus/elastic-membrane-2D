%% interface_comparisons
close all;
clear;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex')


% position = "centre";
% parent_dir = sprintf("/media/michael/newarre/elastic_membrane/%s_tests_8_3_22", position);
parent_dir = "/media/michael/newarre/elastic_membrane/PressureBalanceTests";

% Data directories
dirNames = ["CURVATURE_0-BODYFORCE_0-IFORCE_0", "CURVATURE_1-BODYFORCE_0-IFORCE_0", ...
    "CURVATURE_1-BODYFORCE_1-IFORCE_0", "CURVATURE_1-BODYFORCE_0-IFORCE_1", ...
    "CURVATURE_1-BODYFORCE_1-IFORCE_1"];
dispNames = ["None untransformed", "Curvature transformed", "Curvature and pressure transformed", ...
    "Curvature and surface tension transformed", "All transformed"];

% dirNames = ["CURVATURE_1-BODYFORCE_0-IFORCE_0"];
% dispNames = ["Curvature transformed"];


% Variables
BOX_WIDTH = 6;
xCentre = 3.0;
yCentre = 3.0;
mag = 0.5;


% Exact solution
WExact = @(x) 0.5 * mag * (cos(pi * x / BOX_WIDTH) + 1);

%% Plot convergence
close(figure(5));
figure(5);
hold on;
for k = 1 : length(dirNames)
    % Load in amplitudes file 
    dataDir = sprintf("%s/%s", parent_dir, dirNames(k));
    outputFilename = sprintf("%s/raw_data/amplitudes_8.txt", dataDir);
    amp = readmatrix(outputFilename);
    
    % Plot the covergence 
    ts = amp(:, 1);
    dc = amp(:, 4);
    plot(ts, dc, 'Displayname', dispNames(k), 'Linewidth', 2);
end
xlabel("$t$");
ylabel("change(c, cn)");
grid on;
xlim([0, 5]);
set(gca, 'Yscale', 'log');
legend('Location', 'Southoutside');
set(gcf, 'position', [200 200 1200 900]);

%% Loop over timesteps
for TIMESTEP = 0 : 40
    
    % Loop over cases
    for k = 1 : length(dirNames)
        dataDir = sprintf("%s/%s", parent_dir, dirNames(k));
        
        % Load in solution
        outputFilename = sprintf("%s/raw_data/interface_%d_8.txt", dataDir, TIMESTEP);
        
        % Load in the unsorted matrix A
        A = readmatrix(outputFilename);
        Xs = A(:, 1);
        Ys = A(:, 2);
        
        %% Plot the lab frame
        subplot(1, 2, 1);
        scatter(Xs, Ys - WExact(Xs), 5, 'Displayname', dispNames(k)); 
        hold on;
        
        %% Plot the curvillinear frame
        subplot(1, 2, 2);
        scatter(Xs, Ys, 5, 'Displayname', dispNames(k)); 
        hold on;
    end
    
    %% Plot target solution and membrane in lab frame
    subplot(1, 2, 1);
    % Target solution
    labCurve = @(x, y) (x - xCentre).^2 + (y - yCentre).^2 - 1;
    fimplicit(labCurve, [0 BOX_WIDTH 0 BOX_WIDTH], 'linewidth', 2, 'color', 'black', ...
                'Linestyle', ':', 'Displayname', 'Target solution');
            
    % Membrane
    xs = linspace(0, BOX_WIDTH, 2.^8);
    WVals = WExact(xs);
    scatter(xs, -WVals, 'Displayname', 'Membrane');
%     xline(0, 'color', 'black', 'linewidth', 2);
%     yline(0, 'color', 'black', 'linewidth', 2);
    
    %% Plot target solution and membrane in curvillinear frame
    subplot(1, 2, 2);
    % Target solution
    curvilinearCurve = @(x, y) (x - xCentre).^2 + (y - yCentre - WExact(x)).^2 - 1;
    fimplicit(curvilinearCurve, [0 BOX_WIDTH 0 BOX_WIDTH], 'linewidth', 2, 'color', 'black', ...
                'Linestyle', ':', 'Displayname', 'Target solution');
    
    % Membrane 
    scatter(xs, zeros(size(xs)), 'Displayname', 'Membrane');
%     xline(0, 'color', 'black', 'linewidth', 2);
%     yline(0, 'color', 'black', 'linewidth', 2);


    %% Set lab frame figure properties
    subplot(1, 2, 1);
    hold off;
    xlim([-max(WVals), BOX_WIDTH]);
    ylim([-max(WVals), BOX_WIDTH]);      
    grid on;
    pbaspect([1 1 1]);
    xticks(0:6);
    xlabel('$x$');
    ylabel('$y$');

% legend('Untransformed equations', 'Transformed equations', 'Target solution');
    legend('Location', 'Southoutside');
    title("Lab frame");
    
    %% Set curvillinear frame figure properties
    subplot(1, 2, 2);
    hold off;
    xlim([-max(WVals), BOX_WIDTH]);
    ylim([-max(WVals), BOX_WIDTH]);
    grid on;
    pbaspect([1 1 1]);
    xticks(0:6);
    xlabel('$X$');
    ylabel('$Y$');

%     legend('Untransformed equations', 'Transformed equations', 'Target solution');
    legend('Location', 'Southoutside');
    title("Curvilinear frame");

    
    %% General figure properties
    sgtitle("$t$ =" + num2str(TIMESTEP) , 'Fontsize', 20);
    set(gcf, 'position', [200 200 1200 900]);

    drawnow;
    pause(0.01);
    TIMESTEP

end
