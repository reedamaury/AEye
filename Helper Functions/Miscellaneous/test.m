clear; 
[frame, xx] = obtainframedata(42720, 42734, "Korea_Driving_Data_10.10.10.48_%d.csv");
[xyz,xyz_road] = spherical2cartesian(frame);
plot_xyz(xyz, 1, 100)


[frameboth, xxb] = obtainframedata(1715, 1724, "Both_10.10.10.66_%d.csv");
[xyzb,xyz_roadb] = spherical2cartesian(frameboth);
figure;
plot_xyz(xyzb, 1, 80)
hold on
scatter3(xy_points(:, 1), xy_points(:, 2), zeros(1:length(x), 1)-1.5, "filled", "MarkerFaceColor",'r');
%axis equal;
xlim([0 40]);
ylim([-10 10]);


% [frames_with_noise,frames_with_noise_road] = spherical2cartesian(frame);
% [xyz,xyz_road] = spherical2cartesian(frames_with_noise);
% figure; 
% scatter3(frames_with_noise{1}(:,1),frames_with_noise{1}(:,4),frames_with_noise{1}(:,3),10,frames_with_noise{1}(:,2),'filled');
% plot_frames_with_noise(frames_with_noise_road, 2)
% 
% cov = [1, 1e-1, 1e-6, 0];
% [frames_with_noise_noisy] = add_noise2(frames_with_noise_road, 20, cov, 5);
% plot_frames_with_noise(frames_with_noise_noisy, 2)


