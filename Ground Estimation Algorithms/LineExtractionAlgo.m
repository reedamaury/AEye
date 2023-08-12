% [frameboth, xxb] = obtainframedata(1715, 1724, "Both_10.10.10.66_%d.csv");
% [frameboth, xxb] = obtainframedata(42720, 42734, "Korea_Driving_Data_10.10.10.48_%d.csv");

xyzb = spherical2cartesian(frameboth);

scoor = frameboth{2}(:,1:4);
scoor = scoor(scoor(:,1)<200,:);
points = xyzb{2}(:,1:3); 
frame1 = [scoor points];
frame1(:,2) = [];

% frame1 is now a nx6 matrix, where the columns are range, elevation,
% azimuth, x coord, y coord, z coord, for each point respectively
segments = unique(frame1(:,3)); % lines of constant azimuth
prospects = frame1(frame1(:,2) >= 0, :); % candidate ground points have elevation angle less than 0
figure; plot_xyz_no_i([prospects(:,4),prospects(:,5),prospects(:,6)], 150, "r");
prospects(:,1) = vecnorm([prospects(:,4) prospects(:,5)],2,2); % prospects(:,1).*cosd(prospects(:,2)); % changes range to down range


% Define parameters
num_bins = 400; % number of bins, adjust as needed
max_range = 125; % maximum range, adjust as needed
bin_width = max_range / num_bins;

% Initialize bins
bins = linspace(0, max_range, num_bins+1);

% Initialize cell arrays to hold points for each bin in each segment
segment_bins = cell(length(segments), num_bins);
segment_bins_proto = cell(length(segments), num_bins);
uncertain = zeros(1,6);
mt = 0.25;

% Initialize the labels matrix outside the main loop
labels = zeros(size(prospects, 1), 1);

% Create empty matrices to accumulate ground and nonground points
all_ground_points = [];
all_nonground_points = [];

distance_threshold = 0.03; 
lidar_height = -2.1;

tic;
% Assign points to bins
for i = 1:length(segments)
    segment_points = prospects(prospects(:,3) == segments(i), :);
    prototype_points = [];
    n = 1;
     % Create a labels matrix just for the current segment
    segment_labels = zeros(size(segment_points, 1), 1);

    for j = 1:length(bins)-1
        bin_points = segment_points(segment_points(:,1) >= bins(j) & segment_points(:,1) < bins(j+1), :);
        if ~isempty(bin_points)
            [~, min_idx] = min(bin_points(:,end)); % find point with lowest z-coordinate
            prototype_points(n,:) = bin_points(min_idx,:);
            n = n+1;
        end
    end
    % 
    % m = diff(prototype_points(:,6))./diff(prototype_points(:,1));
    % i = 1; 
    % td =0.15;
    % points_left = size(prototype_points,1)-i+1;
  
    % Assume prototype_points is an nx2 matrix
    % rm = min(prototype_points(:,1));
    % rv = rm - 1;
    % virtual_point = [rv 0 segments(i) rv*cosd(segments(i)) rv*cosd(segments(i)) lidar_height];
    % prototype_points = [virtual_point; prototype_points];
    found = false; % To track if we've found a point within the threshold
    for idx = 1:size(prototype_points, 1)
        z_distance = abs(prototype_points(idx, 6) - lidar_height);
        if z_distance <= 0.5
            found = true;
            prototype_points = prototype_points(idx:end, :);
            break;
        end
    end
    
    if ~found % If no points met the criteria, skip the current iteration
        continue;
    end

    r = size(prototype_points, 1);
    
    if r <= 2 
        break; 
    end
    % Initialize
    lines = [];
    line_points = {};

    ii = 1;
    
    
    while ii <= r - 1
        % Start with 5 points, or fewer if there are not enough left
        num_points = min(n - ii + 1, 4);

        loop_counter = 0; 
        while num_points >= 2
             loop_counter = loop_counter + 1; 

            if loop_counter >= 50 % Check if the counter exceeds the threshold
                warning('Exceeded maximum loop iterations! Breaking out of the loop.');
                break;
            end
            % Fit a line to the points using total least squares regression
            [slope, intercept] = tlsregress(prototype_points(ii:ii+num_points-1, [1,6]));
            
            % Compute the predicted y values
            y_pred = slope * prototype_points(ii:ii+num_points-1, 1) + intercept;
    
            % Compute the residuals
            residuals = prototype_points(ii:ii+num_points-1, 6) - y_pred;
    
            % Compute the RMSE
            rmse = sqrt(mean(residuals.^2));
    
            % Check the conditions && abs(intercept-(-2)) <= 0.6 slope 0.05
            % rmes 0.15
            if abs(slope) <= 0.05 && rmse <= 0.15  && abs(intercept-(-2)) <= 0.6
                % The line meets the conditions, so check the next point
                if ii + num_points <= r && abs(prototype_points(ii+num_points, 6) - (slope * prototype_points(ii+num_points, 1) + intercept)) <= 0.025
                    % The next point is close enough, so add it and continue
                    num_points = num_points + 1 ;
                else
                    % The next point is not close enough, so save the line and move on
                    lines = [lines; slope intercept];
                    line_points = [line_points; prototype_points(ii:ii+num_points-1, 1)'];
                    break
                end
            elseif num_points > 2
                % The line doesn't meet the conditions, so try again with fewer points
                num_points = num_points - 1; % 2
            else 
                % The line still doesn't meet conditions, so move on to
                % next point, but only if slope between current point and
                % it is within a threshold (stops us from moving on to some
                % peak)
               num_points = 0; 
               nn = 1; 
               while ii+nn <= r && abs((prototype_points(ii,6)-prototype_points(ii+nn,6))/(prototype_points(ii,1)-prototype_points(ii+nn,1))) > 0.1
                   nn = nn + 1; 
               end
               ii = ii + nn; 
            end
        end
    
        % Move on to the next point, but only if it is close enough to the
        % last line
       if ~isempty(lines)
            while ii + num_points <= n && abs(prototype_points(ii+num_points, 6) - (lines(end,1) * prototype_points(ii+num_points,1) + lines(end,2))) > 0.05
                ii = ii + 1;
            end
       % else
       %     while i + num_points <= n && abs((prototype_points(ii,6)-prototype_points(ii+1,6))/(prototype_points(ii,1)-prototype_points(ii+1,1))) >= 0.06
       %         ii = ii + 1; 
       %     end
       end

       if ii + num_points > r
           break
       else 
           ii = ii + num_points;
       end
       
    end

    % figure;
    % plot(prototype_points(:,1),prototype_points(:,6), ".b")
    % hold on;
    for mm = 1:size(lines,1)
        % Get the x and y coordinates of the line
        x = line_points{mm};
        y = lines(mm, 1) * x + lines(mm, 2);

        % If this is the first line, extend it backwards to the second prototype point
        if mm == 1
            x_extended = [prototype_points(2,1):0.1:x(1), x];  % generate new x values from the second prototype point to the first x value
            y_extended = lines(mm, 1) * x_extended + lines(mm, 2);
            % Update the x-coordinates for this line
            line_points{mm} = x_extended;
            % plot(x_extended, y_extended, 'g');
        end
        % If this is the last line, extend it forwards to the last prototype point
        if mm == size(lines,1)
            x_extended = [x, x(end):0.1:prototype_points(end,1)];  % generate new x values from the last x value to the last prototype point
            y_extended = lines(mm, 1) * x_extended + lines(mm, 2);
            % Update the x-coordinates for this line
            line_points{mm} = x_extended;
            % plot(x_extended, y_extended, 'g');
        end
        if size(lines,1) ~= 1
            % Plot the line
            % plot(x, y, 'g');
        end

        % If this is not the last line, plot a connecting line to the next line
        if mm < size(lines,1)

            % Get the last point of the current line
            last_point = [x(end), y(end)];

            % Get the first point of the next line
            next_x = line_points{mm+1}(1);
            next_y = lines(mm+1, 1) * next_x + lines(mm+1, 2);
            next_point = [next_x, next_y];

            % Plot the connecting line
            % plot([last_point(1), next_point(1)], [last_point(2), next_point(2)], 'g');
            line_points{mm} = x_extended;
        end

         % Update the segment_labels based on distance to the lines
        for idx = 1:size(segment_points, 1)
            point = segment_points(idx, :);
            d = distance_from_line(point(1), point(6), lines(mm, 1), lines(mm, 2));
            if d < distance_threshold
                segment_labels(idx) = 1; % Label as ground
            end
        end

    end

    % hold off;
     % Extract ground and nonground points based on segment_labels
    ground_points = segment_points(segment_labels == 1, :);
    nonground_points = segment_points(segment_labels == 0, :);

    % Accumulate ground and nonground points for final plotting
    all_ground_points = [all_ground_points; ground_points];
    all_nonground_points = [all_nonground_points; nonground_points];
% % Add a title
% title('Line Extraction on Segment of Constant Azimuth');
% % 
% %% Add axis titles
% xlabel('downrange (m)');
% ylabel('Z (m)');
% zlabel('Z (m)');
% 
end
t=toc
% 
% % After the loop, plot all ground and nonground points on a single figure
figure;
hold on;
plot_xyz_no_i([all_ground_points(:,4), all_ground_points(:,5), all_ground_points(:,6)], 140, "r");
plot_xyz_no_i([all_nonground_points(:,4), all_nonground_points(:,5), all_nonground_points(:,6)], 140, "b");
% Add a legend
legend('Ground Points', 'Non-Ground Points');

% Add a title
title('Line Extraction: Ground vs. Non-Ground Points');

% Add axis titles
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
hold off;
% 
% % Assuming all_ground_points is your matrix with ground points
% 
% Extract ground points with z-coordinate higher than -1.63
new_ground_points = all_ground_points(all_ground_points(:,6) > -1.4, :);
fn= [89.8442 0.431225 -1.72512; 87.5593 -0.11465 -1.528835; 75.7017 -1.94927 -1.5825; 77.6406 -2.27041 -1.62703];
zz = zeros(4,3);
false_negatives = [zz fn];
% Extract points based on the given conditions
points_to_add = prospects(prospects(:,4) > 20 & prospects(:,4) < 20.69 & ...  % x-coordinate condition
                          prospects(:,5) > 3.20692 & prospects(:,5) < 3.93563 & ... % y-coordinate condition
                          prospects(:,6) > -2.13 & prospects(:,6) < -2.16, :); % z-coordinate condition

% false_negatives = [false_negatives; points_to_add];
% % Append these points to the new_ground_points matrix
% new_ground_points = [new_ground_points; false_negatives];
% 
% fp= [62.0401 -3.5505 -1.73577; 62.0435 -3.7675 -1.73622; 62.0236 -3.98371 -1.73604; 43.3567 -3.70626 -1.82482];
% % zz = zeros(4,3);
% false_positives = [zz fp];
% % % Assuming false_positives and false_negatives are your matrices for false positives and false negatives respectively
% 
% % 1. Identify false positives based on z-coordinate condition
% additional_false_positives = all_ground_points(all_ground_points(:,6) > -1.4, :);
% false_positives = [false_positives; additional_false_positives];
% 
% % 2. Plot all points in all_ground_points that aren't false positives or false negatives in red
% figure;
% hold on;
% 
% % Mask to identify points in all_ground_points that are not false positives or false negatives
% not_fp_mask = ~ismember(all_ground_points, false_positives, 'rows');
% not_fn_mask = ~ismember(all_ground_points, false_negatives, 'rows');
% mask = not_fp_mask & not_fn_mask;
% 
% plot_xyz_no_i(all_ground_points(mask, 4:6), 150, "r");
% 
% % 3. Plot all false negatives in green
% plot_xyz_no_i(false_negatives(:, 4:6), 150, "g");
% 
% % 4. Plot all false positives in yellow
% plot_xyz_no_i(false_positives(:, 4:6), 150, "c");
% 
% % Add a legend
% legend('Ground Points', 'False Negatives', 'False Positives');
% 
% % Add a title
% title('Line Extraction: Visualization of Ground Points, False Negatives, and False Positives');
% 
% % Add axis titles
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Z (m)');
% 
% % Display the figure with the above annotations
% hold off;
% 




function d = distance_from_line(x, y, m, c)
    a = -m;
    b = 1;
    d = abs(a*x + b*y - c) / sqrt(a^2 + b^2);
end