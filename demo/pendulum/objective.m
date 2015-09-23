function [obj, objGrad] = objective(t,x,u,p)
% [obj, objGrad] = objective(t,x,u,p)
%
% Computes the objective function (and gradients) for the simple pendulum
%

q = x(1,:);
dq = x(2,:);

if nargout == 1   %Numerical gradients
    obj = autoGen_objective(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
else % Analytic gradients
    [obj,objGrad] = autoGen_objective(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
end

end