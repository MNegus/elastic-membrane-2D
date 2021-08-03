%% cutoff_validation.m
% Compares a case with bending stiffness to one without

% close all
clear;

RUN_TYPE = "decoupled";
parent_directory = sprintf("~/Desktop/cutoff_validation/%s", RUN_TYPE);
MAX_TIMESTEP = 4000;
MAXLEVELS = [8, 9, 10, 11, 12, 13];
legend_entries = strings(length(MAXLEVELS), 1);
for q = 1 : length(MAXLEVELS)
   legend_entries(q) = sprintf("Max level = %d", MAXLEVELS(q)); 
end


%% Coupled-decoupled compare

writerobj = VideoWriter(sprintf("%s_max_level_validation.avi", RUN_TYPE));
writerobj.FrameRate = 10;
open(writerobj);
for k = 0 : 10 : MAX_TIMESTEP
    t = k * 1e-4;
    

    % Loops over two cases
    w_max = 0;
    w_min = 1e5;
    
    w_deriv_max = 0;
    w_deriv_min = 1e5;
    for MAXLEVEL = MAXLEVELS
        directory = sprintf("%s/max_level_%d", parent_directory, MAXLEVEL);
        
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

        subplot(3, 1, 1);
        plot(xs, ws, 'linewidth', 1.5);
        hold on;
        
        % Membrane deriv plot
        membrane_mat = dlmread(sprintf("%s/membrane_outputs/w_deriv_%d.txt", directory, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);
        
        if (max(ws) > w_deriv_max)
            w_deriv_max = max(ws);
        end
        
        if (min(ws) < w_deriv_min)
            w_deriv_min = min(ws);
        end

        subplot(3, 1, 2);
        plot(xs, ws, 'linewidth', 1.5);
        hold on;

        % Pressure plot
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", directory, k));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);

        subplot(3, 1, 3);
        plot(xs, ps, 'linewidth', 1.5);
        hold on;
    end
    w_min
    w_max
    subplot(3, 1, 1);
    hold off;
    xlabel("x");
    ylabel("w(x, t)");
    ylim([min(-1e-10, w_min), max(1e-10, w_max)]);
    legend(legend_entries);
    title("t = " + t);
    
    subplot(3, 1, 2);
    hold off;
    xlabel("x");
    ylabel("w_t(x, t)");
    ylim([min(-1e-10, w_deriv_min), max(1e-10, w_deriv_max)]);
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


