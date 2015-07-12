function x = pwPoly2(tGrid,xGrid,t)
% x = pwPoly2(tGrid,xGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in time
%   t = [1, k] = vector of query times (must be contained within tGrid)
%
% OUTPUTS:
%   x = [m, k] = function value at each query time
%
% NOTES: 
%   If t is out of bounds, then all corresponding values for x are replaced
%   with NaN
%

nGrid = length(tGrid);
if mod(nGrid-1,2)~=0 || nGrid < 3
    error('The number of grid-points must be odd and at least 3');
end

% Figure out sizes
n = floor((length(tGrid)-1)/2);
m = size(xGrid,1);
k = length(t);
x = zeros(m, k);

% Figure out which segment each value of t should be on
tGrid(end) = tGrid(end) + 1e-14; %Hack to include queries on the endpoint
edges = [-inf, tGrid(1:2:end), inf];
[~, bin] = histc(t,edges);

% Loop over each quadratic segment
for i=1:n
    idx = bin==(i+1);
    if sum(idx) > 0
        gridIdx = 2*(i-1) + [1,2,3];
        x(:,idx) = quadInterp(tGrid(gridIdx),xGrid(:,gridIdx),t(idx));
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==1 | bin==(n+2);
x(:,outOfBounds) = nan;

end


function x = quadInterp(tGrid,xGrid,t)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   tGrid = [1, 3] = time grid
%   xGrid = [m, 3] = function grid
%   t = [1, p] = query times, spanned by tGrid
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

% Rescale the query points to be on the domain [-1,1]
t = 2*(t-tGrid(1))/(tGrid(3)-tGrid(1)) - 1; 

% Compute the coefficients:
a = 0.5*(xGrid(:,3) + xGrid(:,1)) - xGrid(:,2);
b = 0.5*(xGrid(:,3)-xGrid(:,1));
c = xGrid(:,2);

% Evaluate the polynomial for each dimension of the function:
p = length(t);
m = size(xGrid,1);
x = zeros(m,p);
tt = t.^2;
for i=1:m
    x(i,:) = a(i)*tt + b(i)*t + c(i);
end

end


