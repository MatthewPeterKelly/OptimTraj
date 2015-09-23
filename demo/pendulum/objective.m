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
    [nState, nTime] = size(x);
    nControl = size(u,1);
    nVars = 1+nState+nControl;
    [obj,Cz] = autoGen_objective(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
    objGrad = zeros(1,nTime,nVars);
    for i=1:nTime
        objGrad(1,i,:) = reshape(Cz(:,i),1,nVars);
    end
end

end