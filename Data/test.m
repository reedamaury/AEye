clear; 
[frame, xx] = obtainframedata(1716, 1724, "Both_10.10.10.66_%d.csv");
frames_with_noise = add_noise1(frame, 1);
 

[xyz,xyz_road] = spherical2cartesian(frames_with_noise);
plot_xyz(xyz_road, 2)
% [frames_with_noise,frames_with_noise_road] = spherical2cartesian(frame);
% figure; 
% scatter3(frames_with_noise{1}(:,1),frames_with_noise{1}(:,4),frames_with_noise{1}(:,3),10,frames_with_noise{1}(:,2),'filled');
% plot_frames_with_noise(frames_with_noise_road, 2)
% 
% cov = [1, 1e-1, 1e-6, 0];
% [frames_with_noise_noisy] = add_noise2(frames_with_noise_road, 20, cov, 5);
% plot_frames_with_noise(frames_with_noise_noisy, 2)

