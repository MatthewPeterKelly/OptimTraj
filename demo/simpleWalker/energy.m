function [E,U,T] = energy(z,p)
% [E,U,T] = energy(z,p)
%
% This function computes the mechanical energy of the simple walker
%
% INPUTS:
%   z = [4,n] = state vector
%   p = parameter struct:
%       .m1 = hip mass
%       .m2 = foot mass
%       .g = gravitational acceleration
%       .l = leg length
%
% OUTPUTS:
%   E = [1,n] = total mechanical energy
%   U = [1,n] = potential energy
%   T = [1,n] = kinetic energy
% 
% NOTES:
%   
%   states:
%       1 = q1 = first link angle
%       2 = q2 = second link angle
%       3 = dq1 = first link angular rate
%       4 = dq2 = second link angular rate
%
%   angles: measured from negative j axis with positive convention
%

q1 = z(1,:);
q2 = z(2,:);
dq1 = z(3,:);
dq2 = z(4,:);

[U,T] = autoGen_energy(q1,q2,dq1,dq2,p.m1,p.m2,p.g,p.l);

E = U+T;

end