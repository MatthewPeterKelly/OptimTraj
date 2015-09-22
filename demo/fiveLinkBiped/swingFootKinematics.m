function [P5, dP5] = swingFootKinematics(q,dq,p)
% [P5, dP5] = swingFootKinematics(q,dq,p)
%
% This function computes the kinematics of the swing foot, for use in the
% step constraints.
%
% 

[P5,dP5] = autoGen_swingFootKinematics(...
    q(1,:), q(2,:), q(4,:), q(5,:),...
    dq(1,:), dq(2,:), dq(4,:), dq(5,:),...
    p.l1, p.l2, p.l4, p.l5);

end