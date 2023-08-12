function [inlier_points, outlier_points, inlierIndices, outlierIndices] = extract_ground_pca(points, distance_threshold, height_threshold)
    % 1. Center the data
    mean_points = mean(points, 1);
    centered_points = points - mean_points;

    % 2. Compute the covariance matrix
    cov_matrix = cov(centered_points);

    % 3. Compute eigenvectors and eigenvalues
    [V, D] = eig(cov_matrix);

    % 4. The eigenvector corresponding to the smallest eigenvalue is the normal of the ground
    [~, min_index] = min(diag(D));
    ground_normal = V(:, min_index);
    
    % Ensure the ground normal is oriented towards the vertical axis
    if ground_normal(3) < 0
        ground_normal = -ground_normal;
    end

    % Use the normal vector to determine inliers and outliers
    distance_to_plane = abs(dot(centered_points, repmat(ground_normal', size(points, 1), 1), 2));
    
    potential_ground_indices = find(distance_to_plane < distance_threshold);
    
    % Apply height-based filtering
    potential_ground_points = points(potential_ground_indices, :);
    heights = potential_ground_points(:, 3) - mean_points(3);
    height_inliers = abs(heights) < height_threshold;
    
    inlierIndices = potential_ground_indices(height_inliers);
    outlierIndices = setdiff(1:size(points, 1), inlierIndices);

    % Extract inliers and outliers
    inlier_points = points(inlierIndices, :);
    outlier_points = points(outlierIndices, :);
end
