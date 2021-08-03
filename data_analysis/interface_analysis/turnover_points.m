function ds = turnover_points(output_range, parent_dir, ...
    plate_tol, bubble_box_height, bubble_box_width)
%TURNOVER_POINTS Reads simulation output to find the turnover points
%   Function to read the interface files, which are of the form
%   "interface_n.txt", where n is an integer in output_range and the files
%   are in data_dir. plate_position indicates the vertical position we
%   expect the plate to be in, which can be re-defined as the zero-vertical
%   component. 

% Matrix to save turnover points in
ds = zeros(length(output_range), 2);

for k = 1 : length(output_range)
    % Name of interface file
    filename ...
        = sprintf('%s/interfaces/interface_%d.txt', ...
            parent_dir, output_range(k))
        
    % Reads in the start and end points of the line segments, with x along
    % the horizontal axis and y along the vertical
    [start_points, end_points] = read_interface_points(filename, true); 
    
    % Finds unique values of y in all the points
    all_points = [start_points; end_points];
    [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
    uniq_points = all_points(uniq_idxs, :);
    
    % Removes all points that are below the plate_tol
    uniq_points = uniq_points(uniq_points(:, 1) > plate_tol, :);
    
    % Removes all points inside the bubble box
    uniq_points = uniq_points((uniq_points(:, 1) > bubble_box_height) ...
        | (uniq_points(:, 2) > bubble_box_width), :);
    
    % Sorts the resulting vector in increasing y order
    [~, sorted_idxs] = sort(uniq_points(:, 1));
    sorted_points = uniq_points(sorted_idxs, :);
    
    % Define x as a continuous function of z using interpolation
    x_interp = @(y) interp1(sorted_points(:, 1), sorted_points(:, 2), y);
    
    % Vertical range to search for turnover point in. Corresponds to the
    % bottom of the droplet up to the middle, assuming that the impact is
    % not far along enough that the turnover point is above the vertical
    % mid-point of the droplet.
    ys = linspace(min(sorted_points(:, 1)), 0.25 * max(sorted_points(:, 1)), 1e5);
    
    % Finds local minima of the r_interp function
    local_mins = islocalmin(x_interp(ys));
    
    if nnz(local_mins) == 0
        % If no local minima are found, set the turnover point to just be
        % at the point with the lowest y
        [turnover_y, idx] = min(sorted_points(:, 1));
        turnover_x = sorted_points(idx, 2);
    else  
        % Else, the turnover point is the local minimum with the highest y 
        % value
        turnover_y = min(ys(local_mins));
        turnover_x = x_interp(turnover_y);
    end
        
    % Saves r and z coordinate of the turnover point
    ds(k, 1) = turnover_x;
    ds(k, 2) = turnover_y;
    
    % Plots the turnover point on a graph
    figure(1);
    plot(ys, x_interp(ys));
    hold on;
    scatter(ds(k, 2), ds(k, 1));
    hold off;
    title(num2str(k));
    
end

end