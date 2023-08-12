% Define the range of angles from 90 to 0 degrees in 0.2 degree intervals
angles = 90:-0.2:0;

% Convert angles to radians for trigonometric functions
angles_rad = deg2rad(angles);

% Define the length of the lines
line_length = 1;

% Define the number of points per line
num_points = 10;

% Create a figure
figure;

% Hold the figure to plot multiple lines
hold on;

% Loop over each angle
for i = 1:length(angles_rad)
    % Calculate the x and y coordinates of each point along the line
    x = linspace(0, line_length * cos(angles_rad(i)), num_points);
    y = linspace(0, line_length * sin(angles_rad(i)), num_points);
    
    % Plot the points along the line
    scatter(x, y, 'filled');
end

% Set the axis equal to maintain the aspect ratio
axis equal;

% Release the figure
hold off;
