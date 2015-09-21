function dz = dynamics(z,u,p)
% dz = dynamics(z,u,p)
%
% Computes the first-order dynamics of the five-link biped, wrapper for
% dynSs.
%
% INPUTS:
%   z = [10,n] = first-order state = [q; dq];
%   u = [5, 1] = input torque vector
%   p = parameter struct

q = z(1:5,:);
dq = z(6:10,:);
ddq = dynSs(q,dq,u,p);
dz = [dq;ddq];

end