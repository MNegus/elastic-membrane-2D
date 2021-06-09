%%
close all
clear;

parent_directory = "/home/michael/scratch/stationary_max_levels";

MAX_LEVELS = [8, 9, 10, 11, 12];
% MAX_LEVELS = 10;
legend_entries = strings(size(MAX_LEVELS));
for q = 1 : length(legend_entries)
    legend_entries(q) = sprintf("Max level = %d", MAX_LEVELS(q));
end

MAX_TIMESTEP = 4000;




%%
for k = 0 : 10 : MAX_TIMESTEP
    t = k * 1e-4;
    for MAX_LEVEL = MAX_LEVELS 
        % Reads in ws
        membrane_mat = dlmread(sprintf("%s/max_level_%d/membrane_outputs/w_%d.txt", parent_directory, MAX_LEVEL, k));
        xs = membrane_mat(:, 1);
        ws = membrane_mat(:, 2);
        subplot(2, 1, 1);
        plot(xs, ws, 'linewidth', 1.5);
        hold on;
        
        pressure_mat = dlmread(sprintf("%s/max_level_%d/pressure_outputs/p_%d.txt", parent_directory, MAX_LEVEL, k));
        xs = pressure_mat(:, 1);
        ps = pressure_mat(:, 2);
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
    
    pause(0.01);
end
