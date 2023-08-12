
[frameboth, xxb] = obtainframedata(1715, 1724, "Both_10.10.10.66_%d.csv");

min_el_angle = -6; %min(frameboth{1}(:,3));

[xyzb,xyz_roadb] = spherical2cartesian(frameboth);

points = xyzb{1}(:,1:3); 
pc = pointCloud(points);
% This function uses the M-estimator SAmple Consensus (MSAC) algorithm to find the plane
% MSAC algorithm is a variant of the RANdom SAmple Consensus (RANSAC) algorithm.
[model, inlierIndices] = pcfitplane(pc, 0.1, [0, 0, 1], 10);  % adjust distanceThreshold as necessary
% abcd(:,fr) = model.Parameters';
 % Option to plot planes 
plane1 = select(pc,inlierIndices);
% figure(fr)
% pcshow(plane1)   

[xyz_120, xyz_points, beam_ellipses] = GroundModel(min_el_angle, 2);

pc_model = pointCloud(xyz_120);

[tform, new_pc] = pcregistericp(plane1, pc_model);

figure;
pcshowpair(pc_model,new_pc, "MarkerSize", 15, "BackgroundColor", "w")
xlim([0 120])
ylim([-8 8])
zlim([-3 1])

% figure;
% plot_xyz(xyzb, 1, 80)
% hold on