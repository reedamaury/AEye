function [xyz_noisy] = add_noise_gaussian(xyz, n, cov, smoothing_window)
% This function adds Gaussian noise to the xyz points, intensity and shot number
% Input 1: xyz - Cell array of xyz points, intensity and shot number
% Input 2: n - number of samples to be added per shot
% Output: xyz_noisy - Cell array of xyz points with added Gaussian noise

    xyz_noisy = cell(size(xyz));
    xcov = cov(1);
    ycov = cov(2);
    zcov = cov(3);
    ivar = cov(4);

    % Iterate through each frame
    for i = 1:length(xyz)
        points = xyz{i};
        shot_numbers = unique(points(:,5));
        
        % Initialize noisy_points
        noisy_points = [];

        % Iterate through each shot
        for j = 1:length(shot_numbers)
            shot_number = shot_numbers(j);
            points_in_shot = points(points(:,5) == shot_number, :);
            
            % Compute means and add Gaussian noise
            mean_xyz = mean(points_in_shot(:, 1:3),1);
            mean_intensity = mean(points_in_shot(:, 4));

            for k = 1:n
                noisy_xyz = mvnrnd(mean_xyz, diag([xcov, ycov, zcov]));
                noisy_intensity = normrnd(mean_intensity, sqrt(ivar));
                noisy_points = [noisy_points; noisy_xyz, noisy_intensity, shot_number];
            end
        end
        
        % Add noisy points to xyz_noisy
        xyz_noisy{i} = [points; noisy_points];
        
        % Smooth the intensity values with a moving average filter
        xyz_noisy{i}(:,4) = smoothdata(xyz_noisy{i}(:,4), 'movmean', smoothing_window);

        disp(['on frame number', num2str(i)])
    end
end
