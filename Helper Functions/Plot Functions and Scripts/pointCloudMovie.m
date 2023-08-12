function pointCloudMovie(xyz, pauseTime)
    % Create a figure
    figure; 

    % Set the colormap
    colormap jet;

    % Loop over all frames
    for i = 1:length(xyz)
        % Get the points for this frame
        points = xyz{i};

        % Create a point cloud object
        ptCloud = pointCloud(points(:, 1:3), 'Color', repmat(points(:, 4), [1 3]));

        % Plot the point cloud
        pcshow(ptCloud, 'MarkerSize', 50);
        title(['Frame ', num2str(i)]);
        colorbar;
        %clim([min(points(:, 4)), max(points(:, 4))]); % Set color limits to min and max intensity

        % Optionally pause
        if nargin > 1
            pause(pauseTime);
        end

        % Optionally, you can enable this line to capture a frame for a movie
        mov(i) = getframe(gcf);
    end

    % Uncomment the following lines to save and play the movie
    VideoWriter(mov, 'PointCloudMovie.avi');
    implay('PointCloudMovie.avi');
end
