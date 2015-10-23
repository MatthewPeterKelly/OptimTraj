function [Fx, Fy] = contactForces(q,dq,ddq,p)
% [Fx, Fy] = contactForces(q,dq,ddq,p)
%
% This function computes the contact forces for the five-link biped
%
% INPUTS:
%   q = [5,n] = configuration
%   dq = [5,n] = rates
%   ddq = [5,n] = accelerations
%   p = parameter struct
%
% OUTPUTS:
%   Fx = [1,n] = horizontal contact force acting on robot
%   Fy = [1,n] = vertical contact force acting on robot
%

[Fx,Fy] = autoGen_contactForce(...
    q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
    dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
    ddq(1,:),ddq(2,:),ddq(3,:),ddq(4,:),ddq(5,:),...
    p.m1, p.m2, p.m3, p.m4, p.m5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5, p.g);

end