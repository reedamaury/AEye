function points_rotated = rotatePlane(points, normal)
% This function rotates a plane with a normal vector of [0, 0, 1] to have
% the same orientation as another plane with a given normal vector.
%
% Inputs:
% - points: an nx3 matrix of [x y z] coordinates representing Plane2
% - normal: a 1x3 vector representing the normal vector of Plane1
%
% Output:
% - points_rotated: an nx3 matrix of [x y z] coordinates representing the
%   rotated Plane2

% Normalize the normal vector
normal = normal / norm(normal);

% Calculate the rotation axis as the cross product of the normal vector and the z-axis
rotation_axis = cross(normal, [0, 0, 1]);

% Calculate the rotation angle as the arccosine of the dot product of the normal vector and the z-axis
rotation_angle = acos(dot(normal, [0, 0, 1]));

% Define the rotation matrix using the axis-angle representation
K = [0 -rotation_axis(3) rotation_axis(2); rotation_axis(3) 0 -rotation_axis(1); -rotation_axis(2) rotation_axis(1) 0];
R = eye(3) + sin(rotation_angle)*K + (1-cos(rotation_angle))*K^2;

% Rotate the points
points_rotated = (R * points')';

end
