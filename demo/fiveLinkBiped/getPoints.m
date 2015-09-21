function [P, G] = getPoints(q,p)
% [P, G] = getPoints(q,p)
%
% This function computes the joint positions (P) and center of mass
% positions (G) for the five-link biped, in configuration given by q with
% physical parameters in p.
%
% INPUTS:
%   q = [5, n] = link configuration
%   p = parameters struct
%
% OUTPUTS:
%   P = [10, n] = joint positions [x;y;x;y;...]
%   G = [10, n] = CoM positions {x;y;x;y;...]
%

q1 = q(1,:);  %stance leg tibia angle
q2 = q(2,:);  %stance leg femur angle
q3 = q(3,:);  %torso angle
q4 = q(4,:);  %swing leg femur angle
q5 = q(5,:);  %swing leg tibia angle

[P,G] = autoGen_getPoints(...
    q1,q2,q3,q4,q5,...
    p.l1 ,p.l2 ,p.l3 ,p.l4 ,p.l5 ,p.c1 ,p.c2 ,p.c3 ,p.c4 ,p.c5);

end