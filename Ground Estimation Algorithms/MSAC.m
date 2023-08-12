% [frameboth, xx] = obtainframedata(1716, 1724, "Both_10.10.10.66_%d.csv");
xyz = spherical2cartesian(frameboth);
data = xyz{1}(:,1:3); 
pc = pointCloud(data);
tic
 % This function uses the M-estimator SAmple Consensus (MSAC) algorithm to find the plane
% % MSAC algorithm is a variant of the RANdom SAmple Consensus (RANSAC) algorithm.
[model, inlierIndices] = pcfitplane(pc, 0.05);  % adjust distanceThreshold as necessary
t = toc;
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
