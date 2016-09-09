function soln = chebyshev(problem)
% soln = chebyshev(problem)
%
% This function transcribes a trajectory optimization problem Chebyshev
% orthogonal polynomials for basis functions. This is an orthogonal
% collocation method, where the entire trajectory is represented as a
% single polynomial. It is for problems where the solution can be
% gaurenteed to be smooth to the same degree as the order of the underlying
% polynomial (nColPts-1).
%
% The technique is described in detail in the paper:
%
%   " A Chebyshev Technique for Solving Nonlinear Optimal Control Problems"
%   ISSS Trans. Automatic Control, 1988
%   by:  Jacques Vlassenbroeck  and  Rene Van Dooren
%
% My implementation for computation of the differentiation matrix,
% quadrature rules, and interpolation are based on the following:
%
%   "Barycentric Lagrange Interpolation"
%   Siam Review, 2004
%   Publisher: Society for Industrial and Applied Mathematics
%   by:  Jean-Paul Berrut  and  Lloyd N. Trefethen
%
%   "Approximation Theory and Approximation Practice"
%   Textbook by Lloyd N. Trefethen
%
%   "Chebfun"  Matlab toolbox
%   Website:  http://www.chebfun.org/
%   by  Lloyd N. Trefethen   et al.
%
% For details on the input and output, see the help file for optimTraj.m
%
% Method specific parameters:
%
%   problem.options.method = 'chebyshev'
%   problem.options.chebyshev = struct with method parameters:
%       .nColPts = number of collocation points
%


%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

nColPts = Opt.chebyshev.nColPts;  %Number of grid points for transcription

% Print out some solver info if desired:
if Opt.verbose > 0
    disp('  -> Transcription via Chebyshev orthogonal collocation');
    fprintf('        nColPts = %d \n', nColPts);
    
end

% Compute the parameters for the ORTHogonal polynomial, in this case the
% Chebyshev polynomial roots, quadrature weights, interpolation weights,
% and the differentiation matrix.
try
    [orth.xx, orth.ww, orth.vv] = chebpts(nColPts);
catch ME
    error('Missing dependency:  chebfun  (http://www.chebfun.org/)  ');
end
orth.D = getDifferentiationMatrix(orth.xx,orth.vv);

% Interpolate the guess at the chebyshev-points for transcription:
guess.tSpan = G.time([1,end]);
guess.time = chebpts(nColPts,guess.tSpan)';
guess.state = interp1(G.time', G.state', guess.time')';
guess.control = interp1(G.time', G.control', guess.time')';

[zGuess, pack] = packDecVar(guess.time, guess.state, guess.control);

% Unpack all bounds:
dummyMatrix = zeros(1,nColPts-2);  %This just needs to be the right size

tLow = [B.initialTime.low, dummyMatrix, B.finalTime.low];
xLow = [B.initialState.low, B.state.low*ones(1,nColPts-2), B.finalState.low];
uLow = B.control.low*ones(1,nColPts);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = [B.initialTime.upp, dummyMatrix, B.finalTime.upp];
xUpp = [B.initialState.upp, B.state.upp*ones(1,nColPts-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nColPts);
zUpp = packDecVar(tUpp,xUpp,uUpp);


%%%% Set up problem for fmincon:

P.objective = @(z)( ...
    myObjective(z, pack, F.pathObj, F.bndObj, orth) );

P.nonlcon = @(z)( ...
    myConstraint(z, pack, F.dynamics, F.pathCst, F.bndCst, orth) );

P.x0 = zGuess;
P.lb = zLow;
P.ub = zUpp;
P.Aineq = []; P.bineq = [];
P.Aeq = []; P.beq = [];
P.options = Opt.nlpOpt;
P.solver = 'fmincon';

%%%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal,exitFlag,output] = fmincon(P);
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack,orth);
nlpTime = toc;

%%%% Store the results:
soln.grid.time = tSoln;
soln.grid.state = xSoln;
soln.grid.control = uSoln;

%%%% Rescale the points:
dSoln = tSoln([1,end]);   %Domain of the final solution
xxSoln = orthScale(orth,dSoln);
soln.interp.state = @(t)( barycentricInterpolate(t', xSoln',xxSoln,orth.vv)' );
soln.interp.control = @(t)( barycentricInterpolate(t', uSoln',xxSoln,orth.vv)' );

soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%% SUB FUNCTIONS %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [z,pack] = packDecVar(t,x,u)
%
% This function collapses the time (t), state (x)
% and control (u) matricies into a single vector
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%
% OUTPUTS:
%   z = column vector of 2 + nTime*(nState+nControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nTime
%       .nState
%       .nControl
%

nTime = length(t);
nState = size(x,1);
nControl = size(u,1);

tSpan = [t(1); t(end)];
xCol = reshape(x, nState*nTime, 1);
uCol = reshape(u, nControl*nTime, 1);

z = [tSpan;xCol;uCol];

pack.nTime = nTime;
pack.nState = nState;
pack.nControl = nControl;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [t,x,u,w] = unPackDecVar(z,pack,orth)
%
% This function unpacks the decision variables for
% trajectory optimization into the time (t),
% state (x), and control (u) matricies
%
% INPUTS:
%   z = column vector of 2 + nTime*(nState+nControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nTime
%       .nState
%       .nControl
%
% OUTPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%   w = [1, nTime] = weights for clenshaw-curtis quadrature
%

nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;
nx = nState*nTime;
nu = nControl*nTime;

[t, w] = orthScale(orth,[z(1),z(2)]);
t = t';
x = reshape(z((2+1):(2+nx)),nState,nTime);
u = reshape(z((2+nx+1):(2+nx+nu)),nControl,nTime);


end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function cost = myObjective(z,pack,pathObj,bndObj,cheb)
%
% This function unpacks the decision variables, sends them to the
% user-defined objective functions, and then returns the final cost
%
% INPUTS:
%   z = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   pathObj = user-defined integral objective function
%   endObj = user-defined end-point objective function
%
% OUTPUTS:
%   cost = scale cost for this set of decision variables
%

[t,x,u,w] = unPackDecVar(z,pack,cheb);

% Compute the cost integral along trajectory
if isempty(pathObj)
    integralCost = 0;
else
    integrand = pathObj(t,x,u);  %Calculate the integrand of the cost function
    integralCost = dot(w,integrand);  %Clenshw-curtise quadrature
end

% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    bndCost = bndObj(t0,x0,tF,xF);
end

cost = bndCost + integralCost;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq] = myConstraint(z,pack,dynFun, pathCst, bndCst, orth)
%
% This function unpacks the decision variables, computes the defects along
% the trajectory, and then evaluates the user-defined constraint functions.
%
% INPUTS:
%   z = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynFun = user-defined dynamics function
%   pathCst = user-defined constraints along the path
%   endCst = user-defined constraints at the boundaries
%
% OUTPUTS:
%   c = inequality constraints to be passed to fmincon
%   ceq = equality constraints to be passed to fmincon
%

[t,x,u] = unPackDecVar(z,pack,orth);


%%%% Enforce the dynamics:

% Analytic differentiation of the trajectory at chebyshev points:
d = t([1,end]);   %Domain of the trajectory
[~,~,D] = orthScale(orth,d);  %Scale the differentiation matrix
dxFun = (D*(x'))';   %Differentiate trajectory

% Derivative, according to the dynamics function:
dxDyn = dynFun(t,x,u);

% Add a constraint that both versions of the derivative must match:
defects = dxFun - dxDyn;


%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,defects, pathCst, bndCst);

end






function [x,w,D] = orthScale(orth,d)
% [x,w,D] = orthScale(orth,d)
%
% This function scales the chebyshev points to an arbitrary interval
%
% INPUTS:
%   xx = chebyshev points on the domain [-1,1]
%   ww = chebysehv weights on the domain [-1,1]
%   d = [low, upp] = new domain
%
% OUTPUTS:
%   x = chebyshev points on the new domain d
%   w = chebyshev weights on the new domain d
%

shift = 0.5*(d(1) + d(2));
scale = 0.5*(d(2) - d(1));

x = scale*orth.xx + shift;

if nargout > 1
    w = orth.ww*scale;
end

if nargout > 2
    D = orth.D/scale;
end

end




function D = getDifferentiationMatrix(x,v,d)
% D = getDifferentiationMatrix(x,v,d)
%
%
% INPUTS:
%   x = [n,1] = vector of roots of the orthogonal polynomial of interest
%   v = [n,1] = vector of barycentric weights corresponding to each root
%   d = [1,2] = domain of the polynomial (optional)
%
% OUTPUTS:
%   D = [n,n] = differentiation matrix such that dy/dx = D*y @ points in x
%
% NOTES:
%   Reference:
%       1) ChebFun  (http://www.chebfun.org/)
%       2) "Barycentric Lagrange Interpolation"   SIAM Review 2004
%           Jean-Paul Berrut and Lloyd N. Trefethen
%
%   Inputs: x and v are typically produced by a call to any of:
%       chebpts, trigpts, legpts, jacpts, lagpts, hermpts, lobpts, radaupts
%

if nargin == 2
    d = [-1,1];
end

n = length(x);
D = zeros(n,n);
for i=1:n
    D(i,:) = (v/v(i))./(x(i)-x);
    D(i,i) = 0;
    D(i,i) = -sum(D(i,:));
end

D = 2*D/(d(2)-d(1));

end





function y = barycentricInterpolate(x,yk,xk,vk)
% y = barycentricInterpolate(x,yk,xk,vk)
%
% Interpolates an orthogonal polynomial using barycentric interpolation
%
% INPUTS:
%   x = [nTime, 1] = vector of points to evaluate polynomial at
%   yk = [nGrid, nOut] = value of the function to be interpolated at each
%       grid point
%   xk = [nGrid, 1] = roots of orthogonal polynomial
%   vk = [nGrid, 1] = barycentric interpolation weights
%
% OUTPUTS:
%   y = [nTime, nOut] = value of the function at the desired points
%
% NOTES:
%   xk and yk should be produced by chebfun (help chebpts)
%

nOut = size(yk,2);
nTime = length(x);
y = zeros(nTime, nOut);

for i=1:nOut
    y(:,i) = bary(x,yk(:,i),xk,vk);
end

end



