function [E,U,T] = acrobotEnergy(z,p)
% [E,U,T] = acrobotEnergy(z,p)
%
% This function computes the mechanical energy of the acrobot.
%
% INPUTS:
%   z = [4,n] = state vector
%   p = parameter struct:
%       .m1 = elbow mass
%       .m2 = wrist mass
%       .g = gravitational acceleration
%       .l1 = length shoulder to elbow
%       .l2 = length elbow to wrist
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

[U,T] = autoGen_acrobotEnergy(q1,q2,dq1,dq2,p.m1,p.m2,p.g,p.l1,p.l2);

E = U+T;

end