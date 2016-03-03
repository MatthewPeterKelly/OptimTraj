function soln = hermiteSimpson(problem)
% soln = hermiteSimpson(problem)
%
% This function transcribes a trajectory optimization problem using the
% Hermite-Simpson (Seperated) method for enforcing the dynamics. It can be
% found in chapter four of Bett's book:
%
%   John T. Betts, 2001
%   Practical Methods for Optimal Control Using Nonlinear Programming
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'hermiteSimpson'
%   problem.options.hermiteSimpson = struct with method parameters:
%       .nGrid = number of grid points to use for transcription
%
% This transcription method is compatable with analytic gradients. To
% enable this option, set:
%   problem.nlpOpt.GradObj = 'on'
%   problem.nlpOpt.GradConstr = 'on'
%
% Then the user-provided functions must provide gradients. The modified
% function templates are as follows:
%
%         [dx, dxGrad] = dynamics(t,x,u)
%                 dx = [nState, nTime] = dx/dt = derivative of state wrt time
%                 dxGrad = [nState, 1+nx+nu, nTime]
%
%         [dObj, dObjGrad] = pathObj(t,x,u)
%                 dObj = [1, nTime] = integrand from the cost function
%                 dObjGrad = [1+nx+nu, nTime]
%
%         [c, ceq, cGrad, ceqGrad] = pathCst(t,x,u)
%                 c = [nCst, nTime] = column vector of inequality constraints  ( c <= 0 )
%                 ceq = [nCstEq, nTime] = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+nu, nTime];
%                 ceqGrad = [nCstEq, 1+nx+nu, nTime];
%
%         [obj, objGrad] = bndObj(t0,x0,tF,xF)
%                 obj = scalar = objective function for boundry points
%                 objGrad = [1+nx+1+nx, 1]
%
%         [c, ceq, cGrad, ceqGrad] = bndCst(t0,x0,tF,xF)
%                 c = [nCst,1] = column vector of inequality constraints  ( c <= 0 )
%                 ceq = [nCstEq,1] = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+1+nx];
%                 ceqGrad = [nCstEq, 1+nx+1+nx];
%

% Each segment needs an additional data point in the middle, thus:
nGrid = 2*problem.options.hermiteSimpson.nSegment+1;

% Print out some solver info if desired:
if problem.options.verbose > 0
    fprintf('  -> Transcription via Hermite-Simpson method, nSegment = %d\n',...
        problem.options.hermiteSimpson.nSegment);
end

%%%% Method-specific details to pass along to solver:

%Simpson quadrature for integration of the cost function:
problem.func.weights = (2/3)*ones(nGrid,1);
problem.func.weights(2:2:end) = 4/3;
problem.func.weights([1,end]) = 1/3;

% Hermite-Simpson calculation of defects:
problem.func.defectCst = @computeDefects;

%%%% The key line - solve the problem by direct collocation:
soln = directCollocation(problem);

% Use method-consistent interpolation
tSoln = soln.grid.time;
xSoln = soln.grid.state;
uSoln = soln.grid.control;
fSoln = problem.func.dynamics(tSoln,xSoln,uSoln);
soln.interp.state = @(t)( pwPoly3(tSoln,xSoln,fSoln,t) );
soln.interp.control = @(t)(pwPoly2(tSoln,uSoln,t));

% Interpolation for checking collocation constraint along trajectory:
%  collocation constraint = (dynamics) - (derivative of state trajectory)
soln.interp.collCst = @(t)( ...
    problem.func.dynamics(t, soln.interp.state(t), soln.interp.control(t))...
    - pwPoly2(tSoln,fSoln,t) );

% Use multi-segment simpson quadrature to estimate the absolute local error
% along the trajectory.
absColErr = @(t)(abs(soln.interp.collCst(t)));
nSegment = problem.options.hermiteSimpson.nSegment;
nState = size(xSoln,1);
quadTol = 1e-12;   %Compute quadrature to this tolerance 
soln.info.error = zeros(nState,nSegment);
for i=1:nSegment
    idx = 2*i + [-1,1];
    soln.info.error(:,i) = rombergQuadrature(absColErr,tSoln([idx(1), idx(2)]),quadTol);
end
soln.info.maxError = max(max(soln.info.error));

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function [defects, defectsGrad] = computeDefects(dt,x,f,dtGrad,xGrad,fGrad)
%
% This function computes the defects that are used to enforce the
% continuous dynamics of the system along the trajectory.
%
% INPUTS:
%   dt = time step (scalar)
%   x = [nState, nTime] = state at each grid-point along the trajectory
%   f = [nState, nTime] = dynamics of the state along the trajectory
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   dtGrad = [2,1] = gradient of time step with respect to [t0; tF]
%   xGrad = [nState,nTime,nDecVar] = gradient of trajectory wrt dec vars
%   fGrad = [nState,nTime,nDecVar] = gradient of dynamics wrt dec vars
%
% OUTPUTS:
%   defects = [nState, nTime-1] = error in dynamics along the trajectory
%   defectsGrad = [nState, nTime-1, nDecVars] = gradient of defects
%

nTime = size(x,2);

iLow = 1:2:(nTime-1);
iMid = iLow + 1;
iUpp = iMid + 1;

xLow = x(:,iLow);
xMid = x(:,iMid);
xUpp = x(:,iUpp);

fLow = f(:,iLow);
fMid = f(:,iMid);
fUpp = f(:,iUpp);

% Mid-point constraint (Hermite)
defectMidpoint = xMid - (xUpp+xLow)/2 - dt*(fLow-fUpp)/4;

% Interval constraint (Simpson)
defectInterval = xUpp - xLow - dt*(fUpp + 4*fMid + fLow)/3;

% Pack up all defects:
defects = [defectMidpoint, defectInterval];


%%%% Gradient Calculations:
if nargout == 2
    
    xLowGrad = xGrad(:,iLow,:);
    xMidGrad = xGrad(:,iMid,:);
    xUppGrad = xGrad(:,iUpp,:);
    
    fLowGrad = fGrad(:,iLow,:);
    fMidGrad = fGrad(:,iMid,:);
    fUppGrad = fGrad(:,iUpp,:);
    
    % Mid-point constraint (Hermite)
    dtGradTerm = zeros(size(xMidGrad));
    dtGradTerm(:,:,1) = -dtGrad(1)*(fLow-fUpp)/4;
    dtGradTerm(:,:,2) = -dtGrad(2)*(fLow-fUpp)/4;
    defectMidpointGrad = xMidGrad - (xUppGrad+xLowGrad)/2 + dtGradTerm + ...
        - dt*(fLowGrad-fUppGrad)/4;
    
    % Interval constraint (Simpson)
    dtGradTerm = zeros(size(xUppGrad));
    dtGradTerm(:,:,1) = -dtGrad(1)*(fUpp + 4*fMid + fLow)/3;
    dtGradTerm(:,:,2) = -dtGrad(2)*(fUpp + 4*fMid + fLow)/3;
    defectIntervalGrad = xUppGrad - xLowGrad + dtGradTerm + ...
        - dt*(fUppGrad + 4*fMidGrad + fLowGrad)/3;
    
    %Pack up the gradients of the defects:
    defectsGrad = cat(2,defectMidpointGrad,defectIntervalGrad);
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functions for interpolation of the control solution
%

function x = pwPoly2(tGrid,xGrid,t)
% x = pwPoly2(tGrid,xGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data,
% given the function value at the edges and midpoint of the interval of
% interest.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in tGrid
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

% Check for any points that are exactly on the upper grid point:
if sum(t==tGrid(end))>0
    x(:,t==tGrid(end)) = xGrid(:,end);
end

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions for interpolation of the state solution
%




function x = pwPoly3(tGrid,xGrid,fGrid,t)
% x = pwPoly3(tGrid,xGrid,fGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data,
% given the function value at the edges and midpoint of the interval of
% interest.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in time
%   fGrid = [m, 2*n-1] = derivative at each grid point in time
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
edges = [-inf, tGrid(1:2:end), inf];
[~, bin] = histc(t,edges);

% Loop over each quadratic segment
for i=1:n
    idx = bin==(i+1);
    if sum(idx) > 0
        kLow = 2*(i-1) + 1;
        kMid = kLow + 1;
        kUpp = kLow + 2;
        h = tGrid(kUpp)-tGrid(kLow);
        xLow = xGrid(:,kLow);
        fLow = fGrid(:,kLow);
        fMid = fGrid(:,kMid);
        fUpp = fGrid(:,kUpp);
        alpha = t(idx) - tGrid(kLow);
        x(:,idx) = cubicInterp(h,xLow, fLow, fMid, fUpp,alpha);
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==1 | bin==(n+2);
x(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
if sum(t==tGrid(end))>0
    x(:,t==tGrid(end)) = xGrid(:,end);
end

end


function x = cubicInterp(h,xLow, fLow, fMid, fUpp,del)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   h = time step (tUpp-tLow)
%   xLow = function value at tLow
%   fLow = derivative at tLow
%   fMid = derivative at tMid
%   fUpp = derivative at tUpp
%   del = query points on domain [0, h]
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

%%% Fix matrix dimensions for vectorized calculations
nx = length(xLow);
nt = length(del);
xLow = xLow*ones(1,nt);
fLow = fLow*ones(1,nt);
fMid = fMid*ones(1,nt);
fUpp = fUpp*ones(1,nt);
del = ones(nx,1)*del;

a = (2.*(fLow - 2.*fMid + fUpp))./(3.*h.^2);
b = -(3.*fLow - 4.*fMid + fUpp)./(2.*h);
c = fLow;
d = xLow;

x = d + del.*(c + del.*(b + del.*a));

end

