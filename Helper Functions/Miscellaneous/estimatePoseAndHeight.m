function [yaw, pitch, roll, lidar_height] = estimatePoseAndHeight(points)
% This function estimates the yaw, pitch, and roll angles and the height of a LiDAR sensor
% given a point cloud. The point cloud is assumed to be an nx3 matrix of [x y z] coordinates.
%
% The yaw, pitch, and roll angles are defined as follows:
% - Yaw is the rotation around the z-axis
% - Pitch is the rotation around the y-axis
% - Roll is the rotation around the x-axis
%
% The height of the LiDAR sensor is estimated as the negative of the minimum z-value in the point cloud.

% Center the point cloud
points_centered = points - mean(points, 1);

% Apply PCA
[coeff,~,~] = pca(points_centered);

% The columns of coeff are the principal components, and form the rotation matrix
rotation_matrix = coeff;

% Find Euler angles
yaw = atan2(rotation_matrix(2,1), rotation_matrix(1,1));
pitch = atan2(-rotation_matrix(3,1), sqrt(rotation_matrix(3,2)^2 + rotation_matrix(3,3)^2));
roll = atan2(rotation_matrix(3,2), rotation_matrix(3,3));

% Convert to degrees
yaw = rad2deg(yaw);
pitch = rad2deg(pitch);
roll = rad2deg(roll);

% Estimate LiDAR height
lidar_height = -min(points(:,3));

end
