% % Define the range of elevation angles from 0 to -6 degrees in 0.15 degree intervals
% elevation_angles = 0:-0.15:-6;
% 
% % Define the range of azimuth angles from -5 to 5 degrees in 0.2 degree intervals
% azimuth_angles = -5:0.2:5;
% 
% % Define the range of distances from 0 to 200
% ranges = linspace(0, 200, length(elevation_angles));
% 
% % Convert angles to radians for trigonometric functions
% elevation_angles_rad = deg2rad(elevation_angles);
% azimuth_angles_rad = deg2rad(azimuth_angles);
% 
% % Create a figure
% figure;
% 
% % Hold the figure to plot multiple ellipses
% hold on;
% 
% % Loop over each combination of elevation angle, azimuth angle, and range
% for i = 1:length(elevation_angles_rad)
%     for j = 1:length(azimuth_angles_rad)
%         % Calculate the minor and major axes of the ellipse
%         b = ranges(i) * tand(0.2/2);
%         a = (ranges(i) * tand(0.1/2)) / sind(elevation_angles(i));
% 
%         % Calculate the x and y coordinates of the center of the ellipse
%         x = ranges(i) * cos(azimuth_angles_rad(j));
%         y = ranges(i) * sin(azimuth_angles_rad(j));
% 
%         % Create an ellipse at the calculated center with the calculated axes
%         ellipse(a, b, 0, x, y);
%     end
% end
% 
% % Set the axis equal to maintain the aspect ratio
% axis equal;
% 
% % Release the figure
% hold off;

% Define the range of elevation angles from 0 to -6 degrees in 0.15 degree intervals
elevation_angles = 0:-0.15:-6;

% Define the range of azimuth angles from -5 to 5 degrees in 0.2 degree intervals
azimuth_angles = -5:0.2:5;

% Convert angles to radians for trigonometric functions
elevation_angles_rad = deg2rad(elevation_angles);
azimuth_angles_rad = deg2rad(azimuth_angles);

% Calculate the range for each elevation angle given that the LiDAR is 1.5m off the ground
ranges = 1.5 ./ tand(-elevation_angles);

% Create a figure
figure;

% Hold the figure to plot multiple ellipses
hold on;
% x = zeros(1, length(elevation_angles_rad));
% y = zeros(1, length(elevation_angles_rad));
% Loop over each combination of elevation angle, azimuth angle, and range
for i = 1:length(elevation_angles_rad)
    for j = 1:length(azimuth_angles_rad)
        % Calculate the minor and major axes of the ellipse
        b = ranges(i) * tand(0.2/2);
        a = (ranges(i) * tand(0.1/2)) / sind(elevation_angles(i));
        
        % Calculate the x and y coordinates of the center of the ellipse
        x = ranges(i) * cos(azimuth_angles_rad(j));
        y = ranges(i) * sin(azimuth_angles_rad(j));

        % Plot the point
        plot(x, y, ".b", "MarkerSize",2);
        
        % Create an ellipse at the calculated center with the calculated axes
        ellipse(a, b, 0, x, y);
    end
end

% Set the axis equal to maintain the aspect ratio
axis equal;

% Release the figure
hold off;


function ellipse(a, b, angle, x0, y0)
    t = linspace(0, 2*pi, 100);
    x = a*cos(t);
    y = b*sin(t);
    r = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    xy = r * [x; y];
    plot(xy(1,:) + x0, xy(2,:) + y0);
end

