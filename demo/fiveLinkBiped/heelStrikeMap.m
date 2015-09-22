function [qNext, dqNext] = heelStrikeMap(q,dq,p)
% [qNext, dqNext] = heelStrikeMap(q,dq,p)
%
% This function computes the heel-strike map and switches the feet. It
% assumes that the robot is right-left symmetric.
%


%%%% Compute the center of mass velocities 
dG = autoGen_comVel(...
    q(1),q(2),q(3),q(4),q(5),...
    dq(1),dq(2),dq(3),dq(4),dq(5),...
    p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5);

%%%% Swap the feet here:
% 1 <--> 5
% 2 <--> 4

q1 = q(5);
q2 = q(4);
q3 = q(3);
q4 = q(2);
q5 = q(1);

dq1m = dq(5);
dq2m = dq(4);
dq3m = dq(3);
dq4m = dq(2);
dq5m = dq(1);

dGx = dG(1:2:end);
dG1mx = dGx(5);
dG2mx = dGx(4);
dG3mx = dGx(3);
dG4mx = dGx(2);
dG5mx = dGx(1);

dGy = dG(2:2:end);
dG1my = dGy(5);
dG2my = dGy(4);
dG3my = dGy(3);
dG4my = dGy(2);
dG5my = dGy(1);

%%%% Dynamics for heel-strike
[MM,ff] = autoGen_dynHs(...
    q1,q2,q3,q4,q5,...
    dq1m,dq2m,dq3m,dq4m,dq5m,...
    dG1mx,dG2mx,dG3mx,dG4mx,dG5mx,...
    dG1my,dG2my,dG3my,dG4my,dG5my,...
    p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5);

%%%%

dqNext = MM\ff;  %Numerically invert the mass matrix
qNext = [q1;q2;q3;q4;q5];  %Pack-up (now flipped) robot configuration

end