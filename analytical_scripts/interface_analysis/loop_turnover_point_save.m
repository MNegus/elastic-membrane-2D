%% loop_turnover_point_save.m
parent_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/gamma_varying/basilisk_data";

ALPHA = 1;

output_range = 1200:4000;

plate_tol = 3e-3;
bubble_box_height = 2 * plate_tol;
bubble_box_width = 2 * plate_tol;

%% Looped
% for BETA = 0
%     for GAMMA = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]
%         param_dir = sprintf("%s/alpha_%g-beta_%g-gamma_%g", parent_dir, ALPHA, BETA, GAMMA);
%         
%         save_turnover_points(output_range, param_dir, plate_tol, bubble_box_height, bubble_box_width);
%     end
% end


%% Stationary case
param_dir = "/media/michael/newarre/elastic_membrane/confirmation_data/stationary_benchmark";
save_turnover_points(output_range, param_dir, plate_tol, bubble_box_height, bubble_box_width);


