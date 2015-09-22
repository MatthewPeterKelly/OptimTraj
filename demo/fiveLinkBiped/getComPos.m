function G = getComPos(q,p)
% G = getComPos(q,p)
%
% This function returns the center of mass position for the entire robot
%
% INPUTS:
%   q = [5,n] = joint configuration
%   p = parameter struct
%
% OUTPUTS:
%   G = [x;y] = robot center of mass location
%

G = autoGen_comPos(...
    q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
    p.m1,p.m2,p.m3,p.m4,p.m5,...
    p.l1,p.l2,p.l3,p.l4,...
    p.c1,p.c2,p.c3,p.c4,p.c5);

end