% Assuming 'points' is a Nx3 matrix from the expression 'points = xyzb{1}(:,1:3);'

% Apply PCA to get initial ground points
z_idx = find(var(points) == min(var(points)));
data_with_index = [points, (1:size(points, 1))'];
mean_z = mean(points(:, z_idx));
std_z = std(points(:, z_idx));
ground_points = data_with_index((data_with_index(:, z_idx) > mean_z - 1.5 * std_z) & (data_with_index(:, z_idx) < mean_z + 1.5 * std_z), :);

% Create a 2D height map from the initial ground points
grid_resolution = 0.5; % e.g., 0.5 meters
max_x = max(ground_points(:,1));
min_x = min(ground_points(:,1));
max_y = max(ground_points(:,2));
min_y = min(ground_points(:,2));

num_cells_x = ceil((max_x - min_x) / grid_resolution);
num_cells_y = ceil((max_y - min_y) / grid_resolution);

height_map = -inf(num_cells_y, num_cells_x); % Initialize with -inf

for i = 1:size(ground_points, 1)
    x_idx = floor((ground_points(i,1) - min_x) / grid_resolution) + 1;
    y_idx = floor((ground_points(i,2) - min_y) / grid_resolution) + 1;
    z_val = ground_points(i,3);
    
    if z_val > height_map(y_idx, x_idx)
        height_map(y_idx, x_idx) = z_val;
    end
end

% Apply morphological opening to the height map
se = strel('disk', 2); % Adjust this based on your data
opened_height_map = imopen(height_map, se);

% Use the difference between the original and opened height map to refine ground classification
diff_height_map = height_map - opened_height_map;
ground_threshold = 0.2; % Adjust this threshold as needed
is_ground = diff_height_map < ground_threshold;

% Extract refined ground and non-ground points from the original point cloud using the mask
refined_ground_points = [];
non_ground_points = [];

for i = 1:size(points, 1)
    x_idx = floor((points(i,1) - min_x) / grid_resolution) + 1;
    y_idx = floor((points(i,2) - min_y) / grid_resolution) + 1;

    % Ensure indices are within valid range
    if x_idx > num_cells_x
        x_idx = num_cells_x;
    elseif x_idx < 1
        x_idx = 1;
    end

    if y_idx > num_cells_y
        y_idx = num_cells_y;
    elseif y_idx < 1
        y_idx = 1;
    end
    
    if is_ground(y_idx, x_idx)
        refined_ground_points = [refined_ground_points; points(i,:)];
    else
        non_ground_points = [non_ground_points; points(i,:)];
    end
end



% Visualization (assuming you have a function called 'plot_points')
figure;
plot_xyz_no_i(refined_ground_points, 120, 'r'); % Plot ground points in red
hold on;
plot_xyz_no_i(non_ground_points, 120, 'b'); % Plot non-ground points in blue
hold off;