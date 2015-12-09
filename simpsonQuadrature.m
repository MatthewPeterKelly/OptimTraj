function x = simpsonQuadrature(fun,tLow,tUpp,n)
%  x = simpsonQuadrature(fun,tLow,tUpp,n)
%
% This function uses simpson quadrature over each of n uniform segments to
% approximate the integral of fun(t) on the interval tLow <= t <= tUpp
%
% INPUTS:
%   fun = function handle
%       f = fun(t)
%           t = [1, nt] = time query points on [a,b]
%           f = [nx, nt] = vector function at query points
%   tLow = scalar lower bound on time
%   tUpp = scalar upper bound on time
%   n = number of segments to divide interval into
%
% OUTPUTS:
%   x = [nx, 1] = integral along each dimension
%

nt = 2*n+1;
t = linspace(tLow,tUpp,nt);
f = fun(t);
nx = size(f,1);
dt = (tUpp-tLow)/n;

% Compute quadrature weights:
w = ones(nx,nt);
w(:,2:2:end) = 4;
w(:,3:2:(end-2)) = 2;

% Compute quadrature:
x = (dt/6)*sum(w.*f,2);

end