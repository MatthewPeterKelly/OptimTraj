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
    [dx,dxGrad] = autoGen_dynamics(t,q,dq,u,p.m,p.g,p.l,p.c,zeros(size(q)));
    [nx,nt] = size(x);
    nu = size(u,1);
    dxGrad = reshape(dxGrad, nx, 1+nx+nu, nt);
end

end