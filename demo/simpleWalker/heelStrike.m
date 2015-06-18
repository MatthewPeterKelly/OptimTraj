function zAfter = heelStrike(zBefore,p)
% zAfter = heelStrike(zBefore,p)
%
% This function computes the step map of the simple walker, including the
% impulsive heel-strike
%
% INPUTS:
%   zBefore = [4,1] = state vector before heel-strike
%   p = parameter struct:
%       .m = leg mass
%       .I = leg moment of inertia about center of mass
%       .d = distance (along leg) from leg CoM to hip
%       .g = gravitational acceleration
%       .l = leg length
%
% OUTPUTS:
%   zAfter = [4,1] = state vector after heel-strike, including foot switch
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

q1 = zBefore(1,:);
q2 = zBefore(2,:);
dq1 = zBefore(3,:);
dq2 = zBefore(4,:);

[q1New,q2New,dq1New,dq2New] = autoGen_heelStrike(q1,q2,dq1,dq2,p.m,p.I,p.d,p.l);

zAfter = [q1New;q2New;dq1New;dq2New];

end