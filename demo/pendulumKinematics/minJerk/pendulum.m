function dv = pendulum(x,v,u)
% dv = pendulum(x,v,u)
%
% Computes the dynamics of a simple pendulum in second-order form
%
% INPUTS: 
%   x = angle
%   v = rate
%   u = torque
% 
% OUTPUTS:
%   dv = acceleration
%

k = 1.0;   % gravity 
b = 0.5;   % damping

dv = -b*v - k*sin(x) + u;

end