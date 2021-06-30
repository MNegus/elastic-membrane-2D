%% gamma_compare.m
% Compares a case with bending stiffness to one without

close all
clear;

parent_directory = "/home/michael/Desktop/cutoff_ranges";
MAX_TIMESTEP = 2000;
XMIN=4;
YMINS=[0, 1, 2, 3];


%% Coupled-decoupled compare

legend_entries = ["ymin = 0", "ymin = 1", "ymin = 2", "ymin = 3"];
writerobj = VideoWriter(sprintf("xmin_%d.avi", XMIN));
open(writerobj);
for k = 1200 : 1 : MAX_TIMESTEP
    t = k * 1e-4;
    

    % Loops over two cases
    w_max = 0;
    w_min = 1e5;
    figure(1);
    for YMIN = YMINS
        directory = sprintf("%s/x_min_%d-y_min_%d", parent_directory, XMIN, YMIN);
        
        % Pressure plot
        pressure_mat = dlmread(sprintf("%s/membrane_outputs/p_%d.txt", directory, k));
        unsorted_xs = pressure_mat(:, 1);
        unsorted_ps = pressure_mat(:, 2);
        [xs, idxs] = sort(unsorted_xs);
        ps = unsorted_ps(idxs);

%         figure(1);
        plot(xs, ps, 'linewidth', 2);
        hold on;
    end
    w_min
    w_max
    
    yline(0, "linestyle", "--", "color", "black", "linewidth", 2);
    hold off;
    
    xlabel("x");
    ylabel("p(x, t)");
    legend([legend_entries, "p = 0"]);
    grid on;
    
    title(sprintf("xmin = %d, t = %.4f", XMIN, t));
    ylim([-2, 14]);
    xlim([0, 1]);
    
    x0=800;
    y0=800;
    width=1000;
    height=600;
    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(gcf);
    writeVideo(writerobj, frame);
    
%     pause(0.01);
end
close(writerobj);

