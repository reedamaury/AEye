% Section 1: Load CSV data and process the data

% Section 1: Load CSV data and process the data
frame=cell(0);
for i=1716:1724
    % Load selected variables from CSV file
    opts = detectImportOptions("Both_10.10.10.66_"+num2str(i)+".csv");
    opts.SelectedVariableNames = ["range_m_","intensity","laserRow","laserCol"];
    X=table2array(readtable("Both_10.10.10.66_"+num2str(i)+".csv",opts));

    % Preprocessing steps
    % 1. Remove entries with range > 200
    X(X(:,1)>200, :) = [];
    % 2. Sort by range
    [~, sortIndex] = sort(X(:, 1));
    X = X(sortIndex, :);
    % 3. Rescale range and intensity
    X(:,1) = X(:,1) / 64;  % Adjust the division factor as needed
    X(:,2) = min(1,X(:,2)/6000); % Adjust the division factor as needed

    % Compute 3D coordinates from range, elevation and azimuth
    r=X(:,1); fee=X(:,3); thet=X(:,4);
    x=r.*cosd(fee).*cosd(thet); y=r.*cosd(fee).*sind(thet); z=r.*sind(fee);
    scatter3(x,y,-z,10,X(:,2),'filled'); colorbar; ylim([-10 10]); zlim([-2 5]); view(-90,35); xlabel("X (m)");ylabel("Y (m)");zlabel("Z (m)")
    % Store the x,y,z coordinates and the intensity in the cell array
    frame{i-1715} = [x, y, z, X(:,2)];
end

% Section 2: Remove everything but the ground from point cloud
% Here, we are assuming the ground points have z value equal or close to zero.
% You may need to adjust this depending on the actual point cloud.
for i=1:length(frame)
    frame{i}(frame{i}(:,3)>0.5,:) = [];
end
% Section 3: Apply Kalman filter to each point cloud
for i=1:length(frame)
    % Load data
    data = frame{i};

    % Parse CSV data into individual variables
    x = data(:,1); y = data(:,2); z = data(:,3);

    % Define grid parameters
    grid_size = 1;  % adjust as needed
    x_grid = floor(min(x)):grid_size:ceil(max(x));
    y_grid = floor(min(y)):grid_size:ceil(max(y));

    % Initialize state estimates and covariances
    % We are estimating two values (z and dz/dt) for each grid cell
    z_estimates = zeros(length(y_grid), length(x_grid), 2);
    z_covariances = zeros(length(y_grid), length(x_grid), 2, 2);

    % Initialize process and measurement noise covariances
    Q = [0.1 0; 0 0.1];  % adjust as needed
    R = 0.1;  % adjust as needed

    % Loop over the data
    for j = 1:length(x)
        % Find which grid cell this point belongs to
        x_cell = find(x_grid <= x(j), 1, 'last');
        y_cell = find(y_grid <= y(j), 1, 'last');

        % Get the current state estimate and covariance for this grid cell
        state_estimate = squeeze(z_estimates(y_cell, x_cell, :));
        P = squeeze(z_covariances(y_cell, x_cell, :, :));

        % Kalman filter prediction step
        % We are using a constant velocity model, so the position prediction is the old position plus the old velocity
        % The velocity prediction is just the old velocity
        state_prediction = [state_estimate(1) + state_estimate(2); state_estimate(2)];
        P = P + Q;

        % Kalman filter update step
        K = P(1, 1) / (P(1, 1) + R);
        state_estimate = state_prediction + K * (z(j) - state_prediction(1));
        P = (1 - K) * P;

        % Save the updated state estimate and covariance
        z_estimates(y_cell, x_cell, :) = state_estimate;
        z_covariances(y_cell, x_cell, :, :) = P;
    end
end
