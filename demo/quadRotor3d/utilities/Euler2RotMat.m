function [R] = Euler2RotMat(theta_x, theta_y, theta_z)
% [R] = Euler2RotMat(theta_x, theta_y, theta_z)
%
% Stand-in for eul2rotm (built-in function MATLAB Robotics System Toolbox)
% Converts Euler angles [rad] to 3x3 Rotation Matrix
%
% This implementation based on:
% S.M.LaValle, Planning Algorithms, Cambridge University Press, 2006, p.98 
% (PDF available here: http://planning.cs.uiuc.edu/ch3.pdf)
% as cited in: http://nghiaho.com/?page_id=846 
%
% Written by Conrad McGreal @ 2019

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