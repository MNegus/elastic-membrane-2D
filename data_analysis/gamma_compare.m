%% bending_compare.m
% Compares a case with bending stiffness to one without

close all
clear;

parent_directory = "~/scratch/gamma_varying";
MAX_TIMESTEP = 4000;
GAMMAS = [0, 0.1, 1, 2];

%% Coupled-decoupled compare

legend_entries = ["Gamma = 0", "Gamma = 0.1", "Gamma = 1", "Gamma = 2"];
writerobj = VideoWriter("gamma_compare.avi");
open(writerobj);
for k = 0 : 20 : MAX_TIMESTEP
    t = k * 1e-4;
    

    % Loops over two cases
    for gamma = GAMMAS
        directory = sprintf("%s/gamma_%g", parent_directory, gamma);
        
        % Membrane plot
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", directory, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);

        subplot(2, 1, 1);
        plot(xs, ws, 'linewidth', 1.5);
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

    
    subplot(2, 1, 1);
    hold off;
    xlabel("x");
    ylabel("w(x, t)");
    legend(legend_entries);
    title("t = " + t);
    
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


