function [p1,p2,dp1,dp2] = kinematics(z,p)
% [p1,p2,dp1,dp2] = kinematics(z,p)
%
% This function computes the mechanical energy of the acrobot.
%
% INPUTS:
%   z = [4,n] = state vector
%   p = parameter struct:
%       .m = leg mass
%       .I = leg moment of inertia about center of mass
%       .d = distance (along leg) from leg CoM to hip
%       .g = gravitational acceleration
%       .l = leg length
%
% OUTPUTS:
%   p1 = [2,n] = position of the hip
%   p2 = [2,n] = position of the swing foot
%   dp1 = [2,n] = velocity of the hip
%   dp2 = [2,n] = velocity of the swing foot
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

[p1,p2,dp1,dp2] = autoGen_kinematics(q1,q2,dq1,dq2,p.d,p.l);

end