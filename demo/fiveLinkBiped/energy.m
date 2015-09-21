function [KE, PE] = energy(q,dq,p)
%[KE, PE] = energy(q,dq,p)
%
% This function computes the mechanical energy for the five-link biped
%
% INPUTS:
%   q = [5,n] = configuration
%   dq = [5,n] = rates
%   p = parameter struct
%
% OUTPUTS:
%   Fx = [1,n] = horizontal contact force acting on robot
%   Fy = [1,n] = vertical contact force acting on robot
%

[KE, PE] = autoGen_energy(...
 q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
    dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
    p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5, p.g);

end