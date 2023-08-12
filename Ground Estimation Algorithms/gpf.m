% Given frame1 data
% Assuming frame1 is a nx6 matrix as mentioned

% Parameters
N_LPR = 45;  % Number of points used to estimate LPR (adjust based on your data)
Th_seeds = 0.3;  % Threshold for points to be considered initial seeds (adjust as needed)
Th_dist = 0.1;  % Threshold distance from the plane (adjust as needed)
N_iter = 10; % Number of iterations (adjust as needed)
N_segs = 10;   % Number of segments in the x-direction (adjust as needed)

% Split the point cloud into N_segs segments based on x-coordinate
x_values = frame1(:,4);  % Extract x coordinates
x_min = min(x_values);
x_max = max(x_values);
segment_width = (x_max - x_min) / N_segs;

all_ground_points = [];

for seg = 1:N_segs
    % Define segment limits
    x_start = x_min + (seg-1) * segment_width;
    x_end = x_start + segment_width;
    
    % Extract points within the current segment
    segment_indices = find(x_values >= x_start & x_values < x_end);
    segment_data = frame1(segment_indices,:);
    
    % Calculate LPR for the segment
    z_values_segment = segment_data(:,6);  % Extract z coordinates for the segment
    sorted_z = sort(z_values_segment);
    LPR = mean(sorted_z(1:N_LPR));

    % Extract initial seed points for the segment
    seed_indices = find(z_values_segment < LPR + Th_seeds);
    seeds = segment_data(seed_indices,:);

    % Iteratively estimate the plane model for the segment
    for i = 1:N_iter
        % Compute covariance matrix for seeds
        mean_seed = mean(seeds(:,4:6));
        covariance_matrix = cov(seeds(:,4:6));
        
        % Extract normal from the covariance matrix
        [V, D] = eig(covariance_matrix);
        [~, idx] = min(diag(D));
        normal = V(:, idx);  % Normal of the plane
        normal = normal(:);
        
        d = -dot(normal, mean_seed);
        
        % For each point in segment, calculate its distance to the plane
        %distances = abs(dot(normal, segment_data(:,4:6), 2) + d) ./ norm(normal);
        % distances = abs(segment_data(:,4:6) * normal' + d) ./ norm(normal);
        % distances = abs(dot(segment_data(:,4:6)', normal) + d) ./ norm(normal);
        distances = abs(segment_data(:,4:6) * normal + d) ./ norm(normal);


        % Update seeds based on the distances
        seed_indices = find(distances < Th_dist);
        seeds = segment_data(seed_indices,:);
    end
    
    % Concatenate the ground points from this segment
    all_ground_points = [all_ground_points; seeds];
end

% All ground points after processing all segments
ground_points = all_ground_points;

% ... (after the prior implementation code)

% Extract non-ground points
all_point_indices = 1:size(frame1,1);
ground_point_indices = ismember(frame1, ground_points, 'rows');
non_ground_indices = ~ground_point_indices;
non_ground_points = frame1(non_ground_indices,:);

% Plot ground points in red
figure;
plot_xyz_no_i(ground_points(:,4:6), 120, 'red');

hold on;

% Plot non-ground points in blue
plot_xyz_no_i(non_ground_points(:,4:6), 120, 'blue');

title('GPF Algorithm: Ground vs Non-Ground Points');
legend('Ground Points', 'Non-Ground Points');
