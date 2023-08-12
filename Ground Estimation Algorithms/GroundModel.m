function [xyz_120, xyz_points, beam_ellipses] = GroundModel(min_el_angle, lidar_height)

    % Define the range of elevation angles from 0 to -6 degrees in 0.15 degree intervals
    elevation_angles = 0:-0.1:min_el_angle;
    
    % Define the range of azimuth angles from -5 to 5 degrees in 0.2 degree intervals
    azimuth_angles = -45:0.2:45;
    
    % Convert angles to radians for trigonometric functions
    elevation_angles_rad = deg2rad(elevation_angles);
    azimuth_angles_rad = deg2rad(azimuth_angles);
    
    % Calculate the range for each elevation angle given that the LiDAR is 1.5m off the ground
    ranges = lidar_height./ sind(-elevation_angles);
    ground_ranges = lidar_height./tand(-elevation_angles);
    
    % Initialize arrays to store ellipse and xy data points
    beam_ellipses = cell(length(elevation_angles_rad), length(azimuth_angles_rad));
    xyz_points = zeros(length(elevation_angles_rad)*length(azimuth_angles_rad), 3);
    % xyz_points = cell(length(elevation_angles_rad), length(azimuth_angles_rad));
    
    % Index for storing xy points
    k = 1;
    
    % Loop over each combination of elevation angle, azimuth angle, and range
    for i = 1:length(elevation_angles_rad)
        for j = 1:length(azimuth_angles_rad)
            % Calculate the minor and major axes of ground ellipse
            b = ranges(i) * tand(0.2/2);
            a = (ranges(i) * tand(0.1/2)) * ((sind(90-0.05))/ sind(elevation_angles(i)+0.05)); % or cosd(0.05)/cosd(angle_of_incidence+0.05)
            
            % Calculate the x and y coordinates of the center of the ellipse
            x = ground_ranges(i) * cos(azimuth_angles_rad(j));
            y = ground_ranges(i) * sin(azimuth_angles_rad(j));
            z = -lidar_height;
            
            % Store the ellipse and xy data points
            beam_ellipses{i, j} = [a, b, 0, x, y];
            xyz_points(k, :) = [x, y, z];
            % xyz_points{i, j} = [x, y, z];
            
            % Increment the index for storing xy points
            k = k + 1;
        end
    end
    
    % % Create a figure to plot the ellipses
    % for i = 1:size(beam_ellipses, 1)
    %     for j = 1:size(beam_ellipses, 2)
    %         % Retrieve the ellipse parameters
    %         a = beam_ellipses{i, j}(1);
    %         b = beam_ellipses{i, j}(2);
    %         angle = beam_ellipses{i, j}(3);
    %         x0 = beam_ellipses{i, j}(4);
    %         y0 = beam_ellipses{i, j}(5);
    % 
    %         % Plot the ellipse
    %         plot_ellipse(a, b, angle, x0, y0);
    %     end
    % end
    % axis equal;
    % hold off;
     
    xyz_120 = xyz_points(xyz_points(:,1)<120 , :);
    
    % % Create a separate figure to plot the xy points
    % figure;
    % l = length(xyz_points(xyz_points(:,1)<100, 2)); 
    % scatter3(xyz_points(xyz_points(:,1)<100 , 1), xyz_points(xyz_points(:,1)<100, 2), zeros(l, 1)-lidar_height, "filled", "MarkerFaceColor",'b');
    % % axis equal;
    % xlim([0 80])
    % ylim([-10 10])


end

% function plot_ellipse(a, b, angle, x0, y0)
%     t = linspace(0, 2*pi, 100);
%     x = a*cos(t);
%     y = b*sin(t);
%     r = [cos(angle) -sin(angle); sin(angle) cos(angle)];
%     xy = r * [x; y];
%     plot(xy(1,:) + x0, xy(2,:) + y0);
% end


