function [c, cGrad] = cst_footVel(x0,xF,p)
% [c, cGrad] = cst_footVel(x0,xF,p)
%
% This function computes the constraint which requires that the swing foot
% leaves the ground at the start of the step and collides with the ground
% at the end of the step.
%
% INPUTS:
%   t0 = time at the start of the trajectory
%   x0 = state at the start of the trajectory
%   tF = time at the end of the trajectory
%   xF = state at the end of the trajectory
%
% OUTPUTS:
%   c = inequality constraint
%   cGrad = gradient of inequality constraints
%

qp = x0(1:5);
qm = xF(1:5);
dqp = x0(6:10);
dqm = xF(6:10);

if nargout == 1   %Numerical gradients
    c = autoGen_cst_footVel(...
        qp(1),qp(2),qp(4),qp(5),...  %Angles "plus" - immediately after heel strike
        qm(1),qm(2),qm(4),qm(5),...  %Angles "minus" - immediately before heel-strike
        dqp(1),dqp(2),dqp(4),dqp(5),...   %Rates "plus" - immediately after heel-strike
        dqm(1),dqm(2),dqm(4),dqm(5),...   %Rates "minus" - immediately before heel-strike
        p.l1, p.l2, p.l4, p.l5);
else  %Analytic gradients
    [c,cGrad] = autoGen_cst_footVel(...
        qp(1),qp(2),qp(4),qp(5),...  %Angles "plus" - immediately after heel strike
        qm(1),qm(2),qm(4),qm(5),...  %Angles "minus" - immediately before heel-strike
        dqp(1),dqp(2),dqp(4),dqp(5),...   %Rates "plus" - immediately after heel-strike
        dqm(1),dqm(2),dqm(4),dqm(5),...   %Rates "minus" - immediately before heel-strike
        p.l1, p.l2, p.l4, p.l5);
end

end