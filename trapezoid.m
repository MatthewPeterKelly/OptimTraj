function soln = trapezoid(problem)
% soln = trapezoid(problem)
%
% This function transcribes a trajectory optimization problem using the
% trapezoid method for enforcing the dynamics. It can be found in chapter
% four of Bett's book:
%
%   John T. Betts, 2001
%   Practical Methods for Optimal Control Using Nonlinear Programming
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'trapezoid'
%   problem.options.trapezoid = struct with method parameters:
%       .nGrid = number of grid points to use for transcription
%
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

% Print out some solver info if desired:
nGrid = problem.options.trapezoid.nGrid;
if problem.options.verbose > 0
    fprintf('  -> Transcription via trapezoid method, nGrid = %d\n',nGrid);
end

%%%% Method-specific details to pass along to solver:

% Quadrature weights for trapezoid integration:
problem.func.weights = ones(nGrid,1);
problem.func.weights([1,end]) = 0.5;

% Trapazoid integration calculation of defects:
problem.func.defectCst = @computeDefects;

%%%% The key line - solve the problem by direct collocation:
soln = directCollocation(problem);

% Use piecewise linear interpolation for the control
tSoln = soln.grid.time;
xSoln = soln.grid.state;
uSoln = soln.grid.control;
soln.interp.control = @(t)( interp1(tSoln',uSoln',t')' );

% Use piecewise quadratic interpolation for the state:
fSoln = problem.func.dynamics(tSoln,xSoln,uSoln);
soln.interp.state = @(t)( bSpline2(tSoln,xSoln,fSoln,t) );

% Interpolation for checking collocation constraint along trajectory:
%  collocation constraint = (dynamics) - (derivative of state trajectory)
soln.interp.collCst = @(t)( ...
    problem.func.dynamics(t, soln.interp.state(t), soln.interp.control(t))...
    - interp1(tSoln',fSoln',t')' );

% Use multi-segment simpson quadrature to estimate the absolute local error
% along the trajectory.
absColErr = @(t)(abs(soln.interp.collCst(t)));
nSegment = nGrid-1;
nState = size(xSoln,1);
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.error = zeros(nState,nSegment);
for i=1:nSegment
    soln.info.error(:,i) = rombergQuadrature(absColErr,tSoln([i,i+1]),quadTol);
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

idxLow = 1:(nTime-1);
idxUpp = 2:nTime;

xLow = x(:,idxLow);
xUpp = x(:,idxUpp);

fLow = f(:,idxLow);
fUpp = f(:,idxUpp);

% This is the key line:  (Trapazoid Rule)
defects = xUpp-xLow - 0.5*dt*(fLow+fUpp);

%%%% Gradient Calculations:
if nargout == 2
    
    xLowGrad = xGrad(:,idxLow,:);
    xUppGrad = xGrad(:,idxUpp,:);
    
    fLowGrad = fGrad(:,idxLow,:);
    fUppGrad = fGrad(:,idxUpp,:);
    
    % Gradient of the defects:  (chain rule!)
    dtGradTerm = zeros(size(xUppGrad));
    dtGradTerm(:,:,1) = -0.5*dtGrad(1)*(fLow+fUpp);
    dtGradTerm(:,:,2) = -0.5*dtGrad(2)*(fLow+fUpp);
    defectsGrad = xUppGrad - xLowGrad + dtGradTerm + ...
        - 0.5*dt*(fLowGrad+fUppGrad);
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = bSpline2(tGrid,xGrid,fGrid,t)
% x = bSpline2(tGrid,xGrid,fGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data.
% The quadratic interpolant is constructed such that the slope matches on
% both sides of each interval, and the function value matches on the lower
% side of the interval.
%
% INPUTS:
%   tGrid = [1, n] = time grid (knot points)
%   xGrid = [m, n] = function at each grid point in tGrid
%   fGrid = [m, n] = derivative at each grid point in tGrid
%   t = [1, k] = vector of query times (must be contained within tGrid)
%
% OUTPUTS:
%   x = [m, k] = function value at each query time
%
% NOTES:
%   If t is out of bounds, then all corresponding values for x are replaced
%   with NaN
%

[m,n] = size(xGrid);
k = length(t);
x = zeros(m, k);

% Figure out which segment each value of t should be on
[~, bin] = histc(t,[-inf,tGrid,inf]);
bin = bin - 1;

% Loop over each quadratic segment
for i=1:(n-1)
    idx = i==bin;
    if sum(idx) > 0
            h = (tGrid(i+1)-tGrid(i));
            xLow = xGrid(:,i);
            fLow = fGrid(:,i);
            fUpp = fGrid(:,i+1);
            delta = t(idx) - tGrid(i);
            x(:,idx) = bSpline2Core(h,delta,xLow,fLow,fUpp);
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==0 | bin==(n+1);
x(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
if sum(t==tGrid(end))>0
    x(:,t==tGrid(end)) = xGrid(:,end);
end

end


function x = bSpline2Core(h,delta,xLow,fLow,fUpp)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   alpha = fraction of the way through the interval
%   xLow = function value at lower bound
%   fLow = derivative at lower bound
%   fUpp = derivative at upper bound
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

%Fix dimensions for matrix operations...
col = ones(size(delta));
row = ones(size(xLow));
delta = row*delta;
xLow = xLow*col;
fLow = fLow*col;
fUpp = fUpp*col;

fDel = (0.5/h)*(fUpp-fLow);
x = delta.*(delta.*fDel + fLow) + xLow;

end



