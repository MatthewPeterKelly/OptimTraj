function soln = chebyshev(problem, defaultOptions)
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
% For details on the input and output, see the help file for trajOpt.m
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
    disp('    ');
end

% Chebyshev points and weights on the default domain
[xx,ww] = chebyshevPoints(nColPts,[-1,1]);
cheb.xx = xx;
cheb.ww = ww;

% Interpolate the guess at the chebyshev-points for transcription:
guess.tSpan = G.time([1,end]);
guess.time = chebyshevPoints(nColPts,guess.tSpan);
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
    myObjective(z, pack, F.pathObj, F.bndObj, cheb) );

P.nonlcon = @(z)( ...
    myConstraint(z, pack, F.dynamics, F.pathCst, F.bndCst, cheb) );

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
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack,cheb);
nlpTime = toc;

%%%% Store the results:

soln.grid.time = tSoln;
soln.grid.state = xSoln;
soln.grid.control = uSoln;

soln.interp.state = @(t)( chebyshevInterpolate(xSoln,t,[tSoln(1),tSoln(end)]) );
soln.interp.control = @(t)( chebyshevInterpolate(uSoln,t,[tSoln(1),tSoln(end)]) );

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

function [t,x,u,w] = unPackDecVar(z,pack,cheb)
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

[t, w] = chebyshevScalePoints(cheb.xx,cheb.ww,[z(1),z(2)]);
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

function [c, ceq] = myConstraint(z,pack,dynFun, pathCst, bndCst, cheb)
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

[t,x,u] = unPackDecVar(z,pack,cheb);


%%%% Enforce the dynamics:

% Analytic differentiation of the trajectory at chebyshev points:
dxFun = chebyshevDerivative(x,t);

% Derivative, according to the dynamics function:
dxDyn = dynFun(t,x,u);

% Add a constraint that both versions of the derivative must match:
defects = dxFun - dxDyn;


%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,defects, pathCst, bndCst);

end




function [x, w] = chebyshevPoints(n,d)
%[x, w] = chebyshevPoints(n,d)
%
% This function is a light-weight version of the function: chebpts.m
% written by Lloyd Trefethen as part of his Chebyshev Polynomial matlab
% toolbox: chebyfun, which can be downloaded from:
% http://www2.maths.ox.ac.uk/chebfun/download/
%
% The algorithm for computing the quadrature weights is also from
% trefethen's toolbox, with the citation:
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
% http://www2.maths.ox.ac.uk/chebfun/and_beyond/programme/slides/wald.pdf
%
% Slight modifications made by Matthew Kelly
% October 27, 2013
% Cornell University
%
% This function returns the n chebyshev points, over the interval d. Error
% checking has not been included on the inputs.
%
% INPUTS:
%   n = [1x1] the desired number of chebyshev points
%   d = [1x2] domain of the polynomial. Default = [-1,1]
%
% OUTPUTS: (2nd-kind chebyshev points and weights)
%   x = [1xn] the n chebyshev points over the interval d
%   w = [1xn] the n chebyshev weights for Clenshaw-Curtis quadrature
%


if n == 1, x = 0; return, end % Special case

%Compute the chebyshev points on the domain [-1,1]:
m = n-1;
x = sin(pi*(-m:2:m)/(2*m));       % Chebyshev points

%Rescale (if necessary):
if nargin~=1
    x = (diff(d)*x + sum(d))/2;
end

%Check if weights are needed:
if nargout==2
    w = chebyshevWeights(n)*diff(d)/2;
end
end


function w = chebyshevWeights(n) % 2nd-kind Chebyshev wieghts
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
% http://www2.maths.ox.ac.uk/chebfun/and_beyond/programme/slides/wald.pdf
if n == 1
    w = 2;
else
    % new
    n = n-1;
    u0 = 1/(n^2-1+mod(n,2));                      % Boundary weights
    L = 0:n-1; r = 2./(1-4*min(L,n-L).^2);        % Auxiliary vectors
    w = [ifft(r-u0) u0];                          % C-C weights
end
end




function Df = chebyshevDerivative(f,x)
%Df = chebyshevDerivative(f,x)
%
% FUNCTION:
%   This function computes the derivative of the chebyshev interpolant
%   at each of the chebyshev nodes.
%
% INPUTS:
%   f = [nState x nPoints] values at each chebyshev node
%   x = chebyshev points
%
% OUTPUTS:
%   Df = the derivative of the chebyshev interpolant at each chebyshev node
%
% NOTES:
%   The derivative at each node is computed by multiplying f by a
%   differentiation matrix. This matrix is [nPoints x nPoints]. If f is a
%   very large order interpolant then computing this matrix may max out the
%   memory available to matlab.
%

n = size(f,2);

%Get the weight vector
w = ones(1,n);
w(2:2:n) = -1;
w([1,end]) = w([1,end])/2;

%First, compute the weighting matrix:
W = (1./w)'*w;

%Next, compute the matrix with inverse of the node distance
X = zeros(n);
for i=1:n
    idx = (i+1):n;
    X(i,idx) = 1./(x(i)-x(idx));
end

%Use the property that this matrix is anti-symetric
X = X - X';

%Compute the i~=j case:
D = W.*X;

%Deal with the i=j case:
D = D - diag(sum(D,2));

%Apply the differentiation matrix
Df = (D*f')';

end



function [x,w] = chebyshevScalePoints(xx,ww,d)
% [x,w] = chebyshevScalePoints(xx,ww,d)
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

x = ((d(2)-d(1))*xx + sum(d))/2;
w = ww*(d(2)-d(1))/2;

end




function [y, Dy, DDy, DDDy] = chebyshevInterpolate(f,t,d)
%[y, Dy, DDy, DDDy] = chebyshevInterpolate(f,t,d)
%
% This function uses barycentric interpolation to evaluate a chebyshev
% polynomial. Technical details are taken from the book: Approximation
% Theory and Approximation Practice by Lloyd Trefethen. The core equation
% that I use can be found as Theorem 5.2 in this book, and is called
% Barycentric Interpolation.
%
% Written by Matthew Kelly
% October 27, 2013
%   Updated: November 8, 2013
% Cornell University
%
% This function computes the value of the chebyshev polynomial that is
% defined by the vector f, at the inputs t. It can also be used to return
% the first and second derivatives of y with respect to t.
%
% INPUTS:
%   f = [KxN] value of the chebyshev polynomial at each of the chebyshev
%       points (these can be computed using chebyshevPoints.m). Each row
%       represents a single set of chebyshev node values for a single
%       state.
%   t = [1xM] vector of inputs to be evaluated (monotonically increasing)
%   d = [1x2] vector specifying the domain of the polynomial
%
% OUTPUTS:
%   y = [KxM] the value of the interpolant at each point in t
%   Dy = first derivative of y with respect to t
%   DDy = second derivative of y with respect to t
%   DDDy = third derivative of y with respect to t
%
% What is happening here, in plain english:
%
%   There are several Chebyshev points (or nodes) that are spread out
%   across the interval d using a special grid spacing. We will call these
%   points x. For each of the points in x, there is a corresponding value
%   of the chebyshev function, we'll call it f.
%
%   Let's say that we want to find the value of the approximation for some
%   input that is not in x. This can be done using a interpolation of the
%   values in f:
%                   y = c(t,x).*f
%
%   Note that the weighting terms are a function of both the desired input
%   and the chebyshev points. It turns out that there is a more stable and
%   efficient way to calculate this interpolation, which is roughly of the
%   form:
%                   y = (k(t,x).*f)/sum(k(t,x))
%
%   This code evaluates k(t,x) which it then uses to evaluate and return y.
%
%
% NOTE - The base algorithm (described above) is not defined for points in
% t that are close to the chebyshev nodes in x. If a point in t matches a
% point in x, then an additional claculation is required. If a point in t
% is very close to a gridpoint, then there is a slight loss of accuracy,
% which is more pronounced.
%

%Check to see if the t vector is on the proper domain:
idxBndFail = t<min(d) | t>max(d);

%Get the chebyshev points
[k,n] = size(f);
x = chebyshevPoints(n,d);
ONE1 = ones(k,1);
ONE2 = ones(1,length(t));

%Loop through each chebyshev node.
num = zeros(k,length(t));
den = zeros(k,length(t));
for i=1:n
    val = ONE1*(1./(t-x(i)));
    if mod(i,2)==1, val=-val; end;
    if i==1 || i==n
        num = num + 0.5*(f(:,i)*ONE2).*val;
        den = den + 0.5*(val);
    else
        num = num + (f(:,i)*ONE2).*val;
        den = den + val;
    end
end

%compute the solution:
y = num./den;

%Check for any values that were too close to nodes and correct them
nanIdx = isnan(y);
if sum(sum(nanIdx))>0
    nanRowIdx = max(nanIdx,[],1);
    y(:,nanRowIdx) = interp1(x',f',t(nanRowIdx)')';
end

%%%% Replace any out-of-bound queries with NaN:
y(:,idxBndFail) = nan;

%%%% Derivative Calculations %%%%
if nargout == 2
    Df = chebyshevDerivative(f,d);
    Dy = chebyshevInterpolate(Df,t,d);
    Dy(:,idxBndFail) = nan;
elseif nargout == 3
    [Df, DDf] = chebyshevDerivative(f,d);
    Dy = chebyshevInterpolate(Df,t,d);
    DDy = chebyshevInterpolate(DDf,t,d);
    Dy(:,idxBndFail) = nan;
    DDy(:,idxBndFail) = nan;
elseif nargout == 4
    [Df, DDf, DDDf] = chebyshevDerivative(f,d);
    Dy = chebyshevInterpolate(Df,t,d);
    DDy = chebyshevInterpolate(DDf,t,d);
    DDDy = chebyshevInterpolate(DDDf,t,d);
    Dy(:,idxBndFail) = nan;
    DDy(:,idxBndFail) = nan;
    DDDy(:,idxBndFail) = nan;
end




end


