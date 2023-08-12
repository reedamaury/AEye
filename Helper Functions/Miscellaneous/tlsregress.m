function [m, b] = tlsregress(points)
    % Calculate the mean of the points
    mean_points = mean(points);
    
    % Subtract the mean from the points
    X = points - mean_points;
    
    % Perform singular value decomposition
    [U, S, V] = svd(X);
    
    % The slope and intercept of the line
    m = -V(1, 2) / V(2, 2);
    b = mean_points(2) - m * mean_points(1);
end
