function [ceq, ceqGrad] = cst_stepLength(xF,p)
% [ceq, ceqGrad] = cst_stepLength(xF,p)
%
% This function computes the equality constraint on the step length (and
% height) for the gait.
%

qm = xF(1:5);

if nargout == 1 %Numerical gradients
    
    ceq = autoGen_cst_steplength(...
        qm(1),qm(2),qm(4),qm(5),...
        p.l1,p.l2,p.l4,p.l5,p.stepLength);
    
else  %Analytic gradients
    [ceq, ceqGrad] = autoGen_cst_steplength(...
        qm(1),qm(2),qm(4),qm(5),...
        p.l1,p.l2,p.l4,p.l5,p.stepLength);
end

end