function dz = dynamics(z,u)
% dz = dynamics(z,u)
%
% This function returns the dynamics of the toy car example. The state of
% the car is its position and orientation in the plane, and the control is
% the rate of change in steering angle.
%
% INPUTS:
%   z = [3, n] = [x;y;θ] = state = [pos; pos; angle];
%   u = [1, n] = steering rate
%
% OUTPUTS:
%  dz = dz/dt
%

% x = z(1,:);   %Dynamics do not depend on x position
% y = z(2,:):   %Dynamics do not depend on y position
th = z(3,:);

dx = cos(th);
dy = sin(th);
dth = u;

dz = [dx;dy;dth];

end