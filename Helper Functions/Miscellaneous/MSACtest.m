% [frameboth, xxb] = obtainframedata(1715, 1724, "Both_10.10.10.66_%d.csv");



[xyzb] = spherical2cartesian(frameboth);


min_el_angle = -6; % min(frameboth{1}(:,3));
% 
points = xyzb{1}(:,1:3); % columns 1-3 are x y z coordinates respectively of each point.
% points4080 = points(points(:,1)>40 & points(:,1)<80, :);
pc = pointCloud(points);
% % This function uses the M-estimator SAmple Consensus (MSAC) algorithm to find the plane
% % MSAC algorithm is a variant of the RANdom SAmple Consensus (RANSAC) algorithm.
tic; 
[model, inlierIndices] = pcfitplane(pc, 0.05);  % adjust distanceThreshold as necessary
t = toc
abcd(1,:) = model.Parameters';

%  % Option to plot planes 
plane1 = select(pc,inlierIndices);
plane_points = plane1.Location;

% Determine the outlier indices
allIndices = 1:size(points, 1);
outlierIndices = setdiff(allIndices, inlierIndices);

% Get the outlier points
outlier_points = points(outlierIndices, :);

% Plotting
figure;
% Plot inliers in red
plot_xyz_no_i(plane_points, 120, "r");
hold on; % This will allow you to plot on the same figure without clearing it

% Plot outliers in blue
plot_xyz_no_i(outlier_points, 120, "b");
hold off; % Done plotting
title('MSAC Algorithm: Ground vs Non-Ground Points');
legend('Ground Points', 'Non-Ground Points');

% 
% 
% [yaw, pitch, roll, lidar_height] = estimatePoseAndHeight(plane_points);
% 
% [xyz_120, xyz_points, beam_ellipses] = GroundModel(min_el_angle, 2.0067);
% 
% xyz4080 = xyz_120(xyz_120(:,1)>40 & xyz_120(:,1)<80, :);
% 
% rotation_point = [0 0 -2.0067];
% 
% xyz_4080_rot = rotatePointCloud(xyz4080, yaw, pitch, roll, rotation_point);
% 
% % normal = abcd(1, 1:3);
% % xyz_120_rot2 = rotatePlane(xyz_120, normal);
% % 
% % rotation_point = [20.0651 0 -1.86519];
% % 
% % xyz_120_rot22 = rotatePointCloud(xyz_120_rot, 2, 0, 0, rotation_point);
% 
% 
% 
% % xyz_120_rot = rotatePointCloud(xyz_120, 0, pitch, roll, rotation_point);
% % 
% % normal = abcd(1, 1:3);
% % xyz_120_rot2 = rotatePlane(xyz_120, normal);
% % 
% % rotation_point = [20.0651 0 -1.86519];
% % 
% % xyz_120_rot22 = rotatePointCloud(xyz_120_rot, 2, 0, 0, rotation_point);
% 
% figure; 
%  plot_xyz_no_i(xyz_4080_rot, 85, "b")
%  hold on 
%  plot_xyz_no_i(plane_points, 85, "r")
% 
