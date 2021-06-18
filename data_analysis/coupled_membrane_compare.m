%%
close all
clear;

parent_directory = "/home/michael/scratch/coupled_max_levels";
MAX_TIMESTEP = 4000;

%% Coupled-decoupled compare

MAX_LEVELS = [12];
legend_entries = strings(2 * length(MAX_LEVELS), 1);
for q = 1 : length(MAX_LEVELS)
    legend_entries(q) = sprintf("Coupled, Max level = %d", MAX_LEVELS(q));
    legend_entries(q + length(MAX_LEVELS)) = sprintf("Decoupled, Max level = %d", MAX_LEVELS(q));
end

writerobj = VideoWriter("coupled_decoupled_compare.avi");
open(writerobj);
for k = 0 : 20 : MAX_TIMESTEP
    t = k * 1e-4;
    
    for COUPLED = [1, 0]
        for MAX_LEVEL = MAX_LEVELS 
            % Membrane plot
            membrane_mat = dlmread(sprintf("%s/max_level_%d-coupled_%d/membrane_outputs/w_%d.txt", parent_directory, MAX_LEVEL, COUPLED, k));
            unsorted_xs = membrane_mat(:, 1);
            unsorted_ws = membrane_mat(:, 2);
            [xs, idxs] = sort(unsorted_xs);
            ws = unsorted_ws(idxs);

            subplot(2, 1, 1);
            plot(xs, ws, 'linewidth', 1.5);
            hold on;

            % Pressure plot
            pressure_mat = dlmread(sprintf("%s/max_level_%d-coupled_%d/membrane_outputs/p_%d.txt", parent_directory, MAX_LEVEL, COUPLED, k));
            unsorted_xs = pressure_mat(:, 1);
            unsorted_ps = pressure_mat(:, 2);
            [xs, idxs] = sort(unsorted_xs);
            ps = unsorted_ps(idxs);

            subplot(2, 1, 2);
            plot(xs, ps, 'linewidth', 1.5);
            hold on;
        end
        
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


%% Coupled-decoupled compare
close all;
MAX_LEVELS = [8, 9, 10, 11, 12];
legend_entries = strings(length(MAX_LEVELS), 1);
for q = 1 : length(MAX_LEVELS)
    legend_entries(q) = sprintf("Coupled, Max level = %d", MAX_LEVELS(q));
end

writerobj = VideoWriter("coupled_maxlevel_compare.avi");
open(writerobj);
for k = 0 : 20 : MAX_TIMESTEP
    t = k * 1e-4;
    
    for MAX_LEVEL = MAX_LEVELS 
        % Membrane plot
        membrane_mat = dlmread(sprintf("%s/max_level_%d-coupled_1/membrane_outputs/w_%d.txt", parent_directory, MAX_LEVEL, k));
        unsorted_xs = membrane_mat(:, 1);
        unsorted_ws = membrane_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ws = unsorted_ws(idxs);

        subplot(2, 1, 1);
        plot(xs, ws, 'linewidth', 1.5);
        hold on;

        % Pressure plot
        pressure_mat = dlmread(sprintf("%s/max_level_%d-coupled_1/membrane_outputs/p_%d.txt", parent_directory, MAX_LEVEL,  k));
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