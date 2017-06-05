function dx = scalarChainIntegrator(x,u)
% dx = scalarChainIntegrator(x,u)
%
% Computes the dynamics for a scalar chain integrator:
%     dx(1) = x(2)
%     dx(2) = x(3)
%     ...
%     dx(n) = u
%

if size(x,1) == 1
    dx = u;
else
    dx = [x(2:end, :); u];
end

end