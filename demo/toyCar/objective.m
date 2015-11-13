function f = objective(z,u)
% f = objective(z,u)
%
% This function computes the integrand of the objective function for the
% toy car problem. The car does not like driving on steep slopes, so it
% will choose the path that minimizes the slope in the direction of travel.
%
% INPUTS:
%   z = [3, n] = [x;y;θ] = state = [pos; pos; angle];
%
% OUTPUTS:
%  f = f(x,y) = integrand of the objective function
%

x = z(1,:);  
y = z(2,:);
th = z(3,:);

% Unit vector in the direction the car is travelling
dx = cos(th);
dy = sin(th);

% Gradient of the height map with respect to x and y
[~, fx, fy] = hills(x,y);

% Slope in the direction of travel:
slope = dx.*fx + dy.*fy;

% Want slopes near zero, so minimize slope squared:
f = slope.^2 + u.^2;

end