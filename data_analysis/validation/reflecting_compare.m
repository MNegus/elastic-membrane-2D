%% coarsen_compare.m
% Compares cases with different levels of coarsening

% close all
clear;

parent_directory = "~/scratch/reflecting_waves";
MAX_TIMESTEP = 4000;
MEMBRANE_RADII = [16, 8, 4];

%% Coupled-decoupled compare

legend_entries = ["Membrane radius = 16", "Membrane radius = 8", "Membrane radius = 4"];
writerobj = VideoWriter("radius_compare.avi");
writerobj.FrameRate = 10;
open(writerobj);
for k = 0 : 10 : MAX_TIMESTEP
    t = k * 1e-4;
    
    w_max = 0;
    w_min = 1e5;
    
    w_deriv_max = 0;
    w_deriv_min = 1e5;
    for RADIUS = MEMBRANE_RADII 
        directory = sprintf("%s/membrane_radius_%d", parent_directory, RADIUS);
        
        % Membrane plot
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_%d.txt", directory, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        if (RADIUS == 16)
            if (max(ws) > w_max)
                w_max = max(ws);
            end

            if (min(ws) < w_min)
                w_min = min(ws);
            end
        end

        subplot(3, 1, 1);
        if (RADIUS == 0)
            plot(xs, ws, 'linewidth', 4, 'color', 0.75 * [1 1 1]);
        else
            plot(xs, ws, 'linewidth', 2);
        end
        hold on;
        
        % Membrane deriv plot
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", directory, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        if (RADIUS == 16)
            if (max(ws) > w_deriv_max)
                w_deriv_max = max(ws);
            end

            if (min(ws) < w_deriv_min)
                w_deriv_min = min(ws);
            end
        end

        subplot(3, 1, 2);
        if (RADIUS == 0)
            plot(xs, ws, 'linewidth', 4, 'color', 0.75 * [1 1 1]);
        else
            plot(xs, ws, 'linewidth', 2);
        end
        hold on;

        % Pressure plot
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", directory, k));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);

        subplot(3, 1, 3);
        if (RADIUS == 0)
            plot(xs, ps, 'linewidth', 4, 'color', 0.75 * [1 1 1]);
        else
            plot(xs, ps, 'linewidth', 2);
        end
        hold on;
    end
    w_min
    w_max
    subplot(3, 1, 1);
    hold off;
    xlabel("x");
    ylabel("w(x, t)");
    ylim([min(-1e-10, 1.2 * w_min), max(1e-10, 1.2 * w_max)]);
    legend(legend_entries);
    title("alpha = 2, beta = 1, gamma = 2. t = " + t);
    
    subplot(3, 1, 2);
    hold off;
    xlabel("x");
    ylabel("w_t(x, t)");
    ylim([min(-1e-10, 1.2 * w_deriv_min), max(1e-10, 1.2 * w_deriv_max)]);
    legend(legend_entries);
    title("t = " + t);
    
    subplot(3, 1, 3);
    hold off;
    xlabel("x");
    ylabel("p(x, t)");
    legend(legend_entries);
    title("t = " + t);
    if (t > 0.125) 
       ylim([-0.5, 2 / t]); 
    end
    x0=400;
    y0=400;
    width=1000;
    height=1000;
    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(gcf);
    writeVideo(writerobj, frame);
    
    pause(0.01);
end
close(writerobj);
