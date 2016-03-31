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

pos = z(1:3,:);
vel = z(4:6,:);
acc = z(7:9,:);
u = u(1:2,:);   %First two dimensions are rotor forces

dz = dynamics([pos;vel], u, p);   %System dynamics
% dPos = dz(1:3,:);
dVel = dz(4:6,:);

ceq = acc - dVel;  % Acceleration-level constraint
c = [];

end