points = xyzb{1}(:,1:3);
% distance_threshold = 0.05; % adjust as necessary
% height_threshold = 0.5; % adjust as necessary, this will further filter points based on height
% 
% [inlier_points, outlier_points, ~, ~] = extract_ground_pca(points, distance_threshold, height_threshold);
% 
% % Plotting
% figure;
% 
% % Plot inliers in red
% plot_xyz_no_i(inlier_points, 120, "r");
% hold on;
% 
% % Plot outliers in blue
% plot_xyz_no_i(outlier_points, 120, "b");
% hold off;

% Assuming 'points' is a Nx3 matrix from the expression 'points = xyzb{1}(:,1:3);'

% Identify the axis corresponding to height (z-axis)
z_idx = find(var(points) == min(var(points)));

% Append an index column to the data for later use
indices = (1:size(points, 1))';
data_with_index = [points, indices];

% Compute the mean and standard deviation for the z-axis
mean_z = mean(points(:, z_idx));
std_z = std(points(:, z_idx));

% Filter points based on height to classify initial ground points
ground_points = data_with_index((data_with_index(:, z_idx) > mean_z - 1.5 * std_z) & (data_with_index(:, z_idx) < mean_z + 1.5 * std_z), :);

% Scale the z-values of the ground points between 0 and 1
min_z = min(ground_points(:, z_idx));
max_z = max(ground_points(:, z_idx));
ground_points(:, z_idx) = (ground_points(:, z_idx) - min_z) / (max_z - min_z);

% Compute the covariance matrix for the ground points
covariance_matrix = cov(ground_points(:, 1:3));

% Perform eigen-decomposition to get eigenvalues and eigenvectors
[eigen_vectors, eigen_values] = eig(covariance_matrix);

% Get the eigenvector corresponding to the smallest eigenvalue
[~, min_eigenvalue_idx] = min(diag(eigen_values));
normal_vector = eigen_vectors(:, min_eigenvalue_idx);

% Project the ground points onto the normal vector to get projection values
projection = ground_points(:, 1:3) * normal_vector;

% Filter ground points based on projection values to further refine the inliers
threshold = 0.2; % As mentioned in the content
inliers = ground_points(abs(projection) < threshold, :);

% Rescale the z-values of the inliers back to their original scale
inliers(:, z_idx) = inliers(:, z_idx) * (max_z - min_z) + min_z;

% Extract the final ground and non-ground points
final_ground_points = inliers(:, 1:3);
final_non_ground_indices = setdiff(data_with_index(:, 4), inliers(:, 4));
final_non_ground_points = points(final_non_ground_indices, :);

% Visualization (assuming you have a function called 'plot_points')
figure;
plot_xyz_no_i(final_ground_points, 120, 'r'); % Plot ground points in red
hold on;
plot_xyz_no_i(final_non_ground_points, 120, 'b'); % Plot non-ground points in blue
hold off;
