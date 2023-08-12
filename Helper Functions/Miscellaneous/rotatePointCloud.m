function points_rotated = rotatePointCloud(points, yaw, pitch, roll, rotation_point)
% This function rotates a point cloud by the specified yaw, pitch, and roll angles.
% The point cloud is assumed to be an nx3 matrix of [x y z] coordinates.
%
% The yaw, pitch, and roll angles are defined as follows:
% - Yaw is the rotation around the z-axis
% - Pitch is the rotation around the y-axis
% - Roll is the rotation around the x-axis
%
% The angles are assumed to be in degrees.
% The rotation is performed about the specified rotation point.

% Convert angles to radians
yaw = deg2rad(yaw);
pitch = deg2rad(pitch);
roll = deg2rad(roll);

% Define the rotation matrices
Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];

% Combine the rotation matrices
R = Rz * Ry * Rx;

% Translate the point cloud so that the rotation point is at the origin
points_translated = points - rotation_point;

% Rotate the points
points_rotated = (R * points_translated')';

% Translate the points back
points_rotated = points_rotated + rotation_point;

end
