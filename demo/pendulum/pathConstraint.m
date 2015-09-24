function [c, ceq, cGrad, ceqGrad] = pathConstraint(t,x,u,p)
% [c, ceq, cGrad, ceqGrad] = pathConstraint(t,x,u,p)
%
% Computes the path constraint
%

q = x(1,:);
dq = x(2,:);

if nargout == 1   %Numerical gradients
     [c,cGrad] = autoGen_pathCst(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
else % Analytic gradients
     [c,cGrad] = autoGen_pathCst(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
end

cGrad = permute(cGrad,[3,1,2]);
ceq = [];
ceqGrad = [];

end