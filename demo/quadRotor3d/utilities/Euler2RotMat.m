function [R] = Euler2RotMat(theta_x, theta_y, theta_z)
% Converts Euler angles [rad] to 3x3 Rotation Matrix

% See: http://nghiaho.com/?page_id=846 

X = [1 0 0
    0 cos(theta_x) -sin(theta_x)
    0 sin(theta_x) cos(theta_x)];

Y = [cos(theta_y) 0 sin(theta_y)
    0 1 0
    -sin(theta_y) 0 cos(theta_y)];

Z = [cos(theta_z) -sin(theta_z) 0
    sin(theta_z) cos(theta_z) 0
    0 0 1];

R = Z*Y*X;