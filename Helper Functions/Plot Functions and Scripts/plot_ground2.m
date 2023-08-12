% Define the range of elevation angles from 0 to -6 degrees in 0.15 degree intervals
elevation_angles = 0:-0.15:-6;

% Define the range of azimuth angles from -5 to 5 degrees in 0.2 degree intervals
azimuth_angles = -30:0.2:30;

% Convert angles to radians for trigonometric functions
elevation_angles_rad = deg2rad(elevation_angles);
azimuth_angles_rad = deg2rad(azimuth_angles);

% Define LiDAR height
lidar_height = 2;

% Calculate the range for each elevation angle given that the LiDAR is 1.5m off the ground
ranges = lidar_height./ sind(-elevation_angles);

% Initialize arrays to store ellipse and xy data points
ellipses = cell(length(elevation_angles_rad), length(azimuth_angles_rad));
xy_points = zeros(length(elevation_angles_rad)*length(azimuth_angles_rad), 2);

% Index for storing xy points
k = 1;

% Loop over each combination of elevation angle, azimuth angle, and range
for i = 1:length(elevation_angles_rad)
    for j = 1:length(azimuth_angles_rad)
        % Calculate the minor and major axes of the ellipse
        b = ranges(i) * tand(0.2/2);
        a = (ranges(i) * tand(0.1/2)) / sind(elevation_angles(i));
        
        % Calculate the x and y coordinates of the center of the ellipse
        x = ranges(i) * cos(azimuth_angles_rad(j));
        y = ranges(i) * sin(azimuth_angles_rad(j));
        
        % Store the ellipse and xy data points
        ellipses{i, j} = [a, b, 0, x, y];
        xy_points(k, :) = [x, y];
        
        % Increment the index for storing xy points
        k = k + 1;
    end
end

% Create a figure to plot the ellipses
figure;
hold on;
for i = 1:size(ellipses, 1)
    for j = 1:size(ellipses, 2)
        if ellipses{i,j}(4) <= 120 
        % Retrieve the ellipse parameters
        a = ellipses{i, j}(1);
        b = ellipses{i, j}(2);
        angle = ellipses{i, j}(3);
        x0 = ellipses{i, j}(4);
        y0 = ellipses{i, j}(5);
        
        % Plot the ellipse
        ellipse(a, b, angle, x0, y0);
        
        end
    end
end
axis equal;
hold on;

% Create a separate figure to plot the xy points
l = length(xy_points(xy_points(:,1)<100, 2)); 
scatter3(xy_points(xy_points(:,1)<100 , 1), xy_points(xy_points(:,1)<100, 2), zeros(l, 1)-0.05, 10, "filled", "MarkerFaceColor",'b');
% axis equal;
xlim([55 80])
ylim([-10 10])
zlim([-1 1.5])


function ellipse(a, b, angle, x0, y0)
    t = linspace(0, 2*pi, 100);
    x = a*cos(t);
    y = b*sin(t);
    r = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    xy = r * [x; y];
    plot(xy(1,:) + x0, xy(2,:) + y0);
end
