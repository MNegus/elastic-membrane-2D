close all;
clear;

% parent_dir = "/media/michael/newarre/elastic_membrane/model_comparison_data/alpha_2-beta_1-gamma_2/dns";
parent_dir = "/media/michael/newarre/elastic_membrane/stationary_membrane"



%% Load in turnover points
turnover_mat = dlmread(sprintf("%s/turnover_points_basilisk.txt", parent_dir));
ts = turnover_mat(:, 1);
ds = turnover_mat(:, 2);
Hs = turnover_mat(:, 3);

timesteps = (1 : length(ts))';
impact_timestep = 1 + sum(ds == 0);

%% Step 2: Plot interface and scatter from ds, Hs
%
jet_areas = zeros(size(timesteps));

plate_tol = 0;
bubble_box_height = 0.05;
bubble_box_width = 0.05;

figure(1);
% for k = 1 : length(output_range)
for timestep_idx = impact_timestep : length(ts)
    
    % Name of interface file
    filename = sprintf('%s/interfaces/interface_%d.txt', parent_dir, ...
        timestep_idx - 1);

    % Jet root position
    d = ds(timestep_idx);
    H = Hs(timestep_idx);

    %% Load in interface points
    % Reads in the start and end points of the line segments, with x along
    % the horizontal axis and y along the vertical
    [start_points, end_points] = read_interface_points(filename, false); 

    % Finds unique values of y in all the points
    all_points = [start_points; end_points];
    [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
    uniq_points = all_points(uniq_idxs, :);

    % Removes all points that are below the plate_tol
    uniq_points = uniq_points(uniq_points(:, 2) > plate_tol, :);

    % Removes all points inside the bubble box
    uniq_points = uniq_points((uniq_points(:, 2) > bubble_box_height) ...
        | (uniq_points(:, 1) > bubble_box_width), :);

    % Removes all points above the jet root
    uniq_points = uniq_points(uniq_points(:, 2) < H, :);

    % Sorts the resulting vector in increasing x order
    [~, sorted_idxs] = sort(uniq_points(:, 1));
    sorted_points = uniq_points(sorted_idxs, :);
    xs = sorted_points(:, 1);
    ys = sorted_points(:, 2);

    %% Integrate interface to get jet area
    if length(xs) > 1
        area = trapz(xs, ys)
        jet_areas(timestep_idx) = area;
    end
    
    %% Plots the interface 
%     figure(1);
    plot(sorted_points(:, 1), sorted_points(:, 2));
    hold on;
    scatter(ds(timestep_idx), Hs(timestep_idx));
    hold off;
    title(num2str(timestep_idx));
    %     ylim([0, 0.05]);
    drawnow;

end

%% Plot the jet area
figure(2);

plot(ts, jet_areas);

%% Save the jet area
save_mat = zeros(length(timesteps), 2);
save_mat(:, 1) = ts;
save_mat(:, 2) = jet_areas;
save_mat
writematrix(save_mat, sprintf("%s/jet_area.txt", parent_dir));