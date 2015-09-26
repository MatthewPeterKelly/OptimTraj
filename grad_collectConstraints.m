function [c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,defects, defectsGrad, pathCst, bndCst, gradInfo)
% [c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,defects, defectsGrad, pathCst, bndCst, gradInfo)
%
% TrajOpt utility function.
%
% Collects the defects, calls user-defined constraints, and then packs
% everything up into a form that is good for fmincon. Additionally, it
% reshapes and packs up the gradients of these constraints.
%
% INPUTS:
%   t = time vector
%   x = state matrix
%   u = control matrix
%   defects = defects matrix
%   pathCst = user-defined path constraint function
%   bndCst = user-defined boundary constraint function
%
% OUTPUTS:
%   c = inequality constraint for fmincon
%   ceq = equality constraint for fmincon
%

ceq_dyn = reshape(defects,numel(defects),1);
ceq_dynGrad = grad_flattenPathCst(defectsGrad);

%%%% Compute the user-defined constraints:
if isempty(pathCst)
    c_path = [];
    ceq_path = [];
    c_pathGrad = [];
    ceq_pathGrad = [];
else
    [c_path, ceq_path, c_pathGradRaw, ceq_pathGradRaw] = pathCst(t,x,u);
    c_pathGrad = grad_flattenPathCst(grad_reshapeContinuous(c_pathGradRaw,gradInfo));
    ceq_pathGrad = grad_flattenPathCst(grad_reshapeContinuous(ceq_pathGradRaw,gradInfo));
end
if isempty(bndCst)
    c_bnd = [];
    ceq_bnd = [];
    c_bndGrad = [];
    ceq_bndGrad = [];
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    [c_bnd, ceq_bnd, c_bndGradRaw, ceq_bndGradRaw] = bndCst(t0,x0,tF,xF);
    c_bndGrad = grad_reshapeBoundary(c_bndGradRaw,gradInfo);
    ceq_bndGrad = grad_reshapeBoundary(ceq_bndGradRaw,gradInfo);
end

%%%% Pack everything up:
c = [c_path;c_bnd];
ceq = [ceq_dyn; ceq_path; ceq_bnd];

cGrad = [c_pathGrad;c_bndGrad]';
ceqGrad = [ceq_dynGrad; ceq_pathGrad; ceq_bndGrad]';


end



%%%% Sub-Functions %%%%



function C = grad_flattenPathCst(CC)
%
% This function takes a path constraint and reshapes the first two
% dimensions so that it can be passed to fmincon
%
if isempty(CC)
    C = [];
else
    [n1,n2,n3] = size(CC);
    C = reshape(CC,n1*n2,n3);
end

end