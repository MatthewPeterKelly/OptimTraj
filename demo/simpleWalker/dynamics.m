function dz = dynamics(z,u,p)
% dz = dynamics(z,u,p)
%
% This function computes the dynamics of the acrobot: double pendulum, two
% point masses, torque motor between links, no friction.
%
% INPUTS:
%   z = [4,n] = state vector
%   u = [1,n] = hip torque
%   p = parameter struct:
%       .m = leg mass
%       .I = leg moment of inertia about center of mass
%       .d = distance (along leg) from leg CoM to hip
%       .g = gravitational acceleration
%       .l = leg length
%
% OUTPUTS:
%   dz = [4,1] = dz/dt = time derivative of the state
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

[ddq1,ddq2] = autoGen_dynamics(q1,q2,dq1,dq2,u,p.d,p.m,p.I,p.g,p.l);

dz = [dq1;dq2;ddq1;ddq2];

end