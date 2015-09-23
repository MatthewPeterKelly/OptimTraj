function [dx, dxGrad] = dynamics(t,x,u,p)
% [dx, dxGrad] = dynamics(t,x,u,p)
%
% Computes the dynamics (and gradients) for the simple pendulum
%

q = x(1,:);
dq = x(2,:);

if nargout == 1   %Numerical gradients
    dx = autoGen_dynamics(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
else % Analytic gradients
    [dx,Fz] = autoGen_dynamics(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
    
    [nState, nTime] = size(x);
    nControl = size(u,1);
    dxGrad = zeros(nState,nState+nControl,nTime);
    for i=1:nTime
        dxGrad(:,:,i) = Fz(:,i);
    end
end

end