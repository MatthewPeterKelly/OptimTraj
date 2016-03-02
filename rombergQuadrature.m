function [x, err] = rombergQuadrature(fun,tSpan,tol)
% [x, err] = rombergQuadrature(fun,tSpan,tol)
%
% Compute the integral(fun), over the domain tSpan, to an accuracy of tol,
% using Romberg quadrature. Fully vectorized.
%
% Good for high-accuracy quadrature over smooth vector functions.
%
% If necessary, this function will automatically sub-divide the interval to
% achieve the desired accuracy. This should only occur when fun is stiff or
% non-smooth.
%
% INPUTS:
%   fun = vector function to be integrated
%       dx = fun(t)
%           t = [1, nt] = time vector
%           dx = [nx, nt] = function value at each point in t
%   tSpan = [tLow, tUpp] = time span (domain) for integration
%   tol = [nx,1] = desired error tolerance along each dimension
%
% OUTPUT:
%   x = [nx,1] = integral along each dimension
%   err = [nx, 1] = error estimate along each dimension
%
% NOTES:
%   algorithm from:
%   http://www.math.usm.edu/lambers/mat460/fall09/lecture29.pdf
%

a = tSpan(1);
b = tSpan(2);
h = b-a;
f0 = fun(tSpan);
nx = size(f0,1);

if isscalar(tol)
    tol = tol*ones(nx,1);
end

nMin = 4;   % Complete at least this many iterations
nMax = 12;  % Maximum iteration count
nSubSegment = 4;  % Sub-divide segment if max iteration is reached

T = zeros(nMax,nMax,nx);  % (iteration, extrapolation, dimension)
Err = zeros(nMax,nx); %(extrapolation, dimension)

T(1,1,:) = 0.5*h*sum(f0,2);
h = h/2;

converged = false;
for j=2:nMax
    
    % Trapazoid method
    i = 1:2^(j-2);
    t = a+h*(2*i-1);
    f = fun(t);
    T(j,1,:) = 0.5*T(j-1,1,:) + reshape(h*sum(f,2),1,1,nx);
    
    % Richardson extrapolation
    for k=2:j
        T(j,k,:) = T(j,k-1,:) + ...
            ( 1/(4^(k-1)-1) )*( T(j,k-1,:) - T(j-1,k-1,:) );
    end
    h = h/2;
    
    % Check the error estimate for convergence
    Err(j,:) = T(j,k,:) - T(j-1,k-1,:);
    err = abs(Err(j,:))';
    if all( err < tol ) && j > nMin
        converged = true;
        break;
    end
end

x = reshape(T(j,j,:),nx,1);

if ~converged  %Then non-smooth problem. Recursively sub-divide interval.
    time = linspace(tSpan(1), tSpan(2), nSubSegment+1);
    x = zeros(nx,1);
    err = zeros(nx,1);
    for i=1:nSubSegment
        [xTmp, errTmp] = rombergQuadrature(fun,time([i, i+1]),tol);
        x = x + xTmp;
        err = max(err,errTmp);
    end
end

end
