function dv = pendulum(x,v,u,param)
% dv = pendulum(x,v,u,param)
%
% Computes the dynamics of a simple pendulum in second-order form
%
% INPUTS: 
%   x = angle
%   v = rate
%   u = torque
%   param = struct with physical parameters
%       .k = gravity torque constant
%       .b = viscous damping constant
% OUTPUTS:
%   dv = acceleration
%

k = param.k;   % gravity 
b = param.b;   % damping

dv = -b*v - k*sin(x) + u;

end