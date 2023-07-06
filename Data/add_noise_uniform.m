function [frames_with_noise] = add_noise_uniform(frames, n)
    % Create a new cell array to hold frames with noise
    frames_with_noise = cell(size(frames));
    
    % Concatenate all frames to calculate global range_max, intensity_mean, and intensity_std
    all_data = vertcat(frames{:});
    range_max = max(all_data(:, 1));
    % intensity_mean = mean(all_data(:, 2));
    % intensity_std = std(all_data(:, 2));
    
    % Iterate over each frame
    for i = 1:length(frames)
        % Get current frame
        frame = frames{i};
        
        % Find unique shot numbers
        shot_numbers = unique(frame(:,5));

        % Calculate mean and standard deviation of intensity for noise generation
        intensity_mean = mean(frame(:,2));
        intensity_std = mean(frame(:,2));
        
        % Initialize an empty matrix to hold the frame with noise
        frame_with_noise = [];
        
        % Iterate over each shot number
        for shot_number = shot_numbers'
            % Get the data for the current shot number
            shot_data = frame(frame(:,5) == shot_number, :);
            
            % Generate 'n' samples of noise for range and intensity
            noise_range = range_max * rand(n, 1);
            noise_intensity = normrnd(intensity_mean, intensity_std, [n, 1]);
            
            % Create noise data with the same elevation, azimuth, and shot number as the original data
            noise_data = [noise_range, noise_intensity, ones(n,1)*shot_data(1,3:5)];
            
            % Add the original shot data and noise data to the frame with noise
            frame_with_noise = [frame_with_noise; shot_data; noise_data];
        end
        
        % Add the frame with noise to the frames_with_noise cell array
        frames_with_noise{i} = frame_with_noise;
    end
end
