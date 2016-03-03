function [c, ceq] = dynPathCst(z,u,p)
% [c, ceq] = dynPathCst(z,u,p)
%
% This function implements the quad-rotor dynamics as a path constraint.
%
% INPUTS:
%   z = [6, n] = [x; y; q; dx; dy; dq] = state matrix
%   u = [2, n] = [u1; u2] = control matrix
%   p = parameter struct:
%       .g = acceleration due to gravity
%       .d = half distance between rotors
%       .m = mass of each rotor (half-mass of the quad-rotor)
%
% OUTPUTS:
%   c = []
%   ceq = [3,n] = equality constraint on acceleration due to motors
%

dz = dynamics(z(1:6,:), u(1:2,:), p);   %System dynamics
ceq = dz(4:6,:) - z(7:9,:);  % Acceleration-level constraint
c = [];

end