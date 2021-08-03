%% gamma_compare.m
% Compares a case with bending stiffness to one without

% close all
clear;

parent_directory = "/home/michael/scratch/initial_membrane_tests/beta_0_gamma_varying";
MAX_TIMESTEP = 4000;
GAMMAS = [0.01, 0.1, 1, 10, 100];

%% Coupled-decoupled compare

legend_entries = ["gamma = 0.01", "gamma = 0.1", "gamma = 1", "gamma = 10", "gamma = 100"];
writerobj = VideoWriter("gamma_compare_beta_0.avi");
open(writerobj);
for k = 0 : 20 : MAX_TIMESTEP
    t = k * 1e-4;
    

    % Loops over two cases
    w_max = 0;
    w_min = 1e5;
    for gamma = GAMMAS
        directory = sprintf("%s/gamma_%g", parent_directory, gamma);
        
        % Membrane plot
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", directory, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        if (max(ws) > w_max)
            w_max = max(ws);
        end
        
        if (min(ws) < w_min)
            w_min = min(ws);
        end

        subplot(2, 1, 1);
        plot(xs, 4 * ws, 'linewidth', 1.5);
        hold on;

        % Pressure plot
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", directory, k));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);

        subplot(2, 1, 2);
        plot(xs, ps, 'linewidth', 1.5);
        hold on;
    end
    w_min
    w_max
    subplot(2, 1, 1);
    hold off;
    xlabel("x");
    ylabel("w(x, t)");
    ylim([min(-1e-10, 4 * w_min), max(1e-10, 4 * w_max)]);
    legend(legend_entries);
    title("alpha = 2, gamma = 2, t = " + t);
    
    subplot(2, 1, 2);
    hold off;
    xlabel("x");
    ylabel("p(x, t)");
    legend(legend_entries);
    title("t = " + t);
    
    x0=800;
    y0=800;
    width=1000;
    height=1000;
    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(gcf);
    writeVideo(writerobj, frame);
    
    pause(0.01);
end
close(writerobj);

