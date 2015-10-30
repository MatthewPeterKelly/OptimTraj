function soln = multiCheb(problem)
% soln = multiCheb(problem)
%
% DEPRICATED
%
%
% *************************************************************************
% This file is no longer used, and is preserved for reference only. The
% numerical methods for connecting segments are not the most stable,
% particularily for low-order polynomials. This file will later be replaced
% with HP orthogonal collocation, based on Legendre polynomials.
% *************************************************************************
%
%
% This function transcribes a trajectory optimization problem Chebyshev
% orthogonal polynomials for basis functions. This is an orthogonal
% collocation method. This method is similiar to Chebyshev, except that
% here I break the trajectory into several segments, rahter than just one.
%
% The technique is similar to the one described in detail in the paper:
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
%   problem.options.method = 'multiCheb'
%   problem.options.multiCheb = struct with method parameters:
%       .nColPts = number of collocation points in each trajectory segment
%       .nSegment = number of segments to break the trajectory into
%
%
% *************************************************************************
%       DEPRICATED
% *************************************************************************
%


%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

nColPts = Opt.multiCheb.nColPts;  %Number of collocation points in each segment
nSegment = Opt.multiCheb.nSegment;  %Fraction of the duration spent in each segment

% Print out some solver info if desired:
if Opt.verbose > 0
    disp('  -> Transcription via Multiple-segment Chebyshev orthogonal collocation');
    disp('    ');
end

% This method seems to fail if a low-order polynomial is used.
% It gives reasonable solutions for medium-high order polynomials
if nColPts < 6
    disp('    WARNING: using fewer than six collocation points per interval can lead to numerical problems!');
end

% Chebyshev points and weights on the default domain
[xx,ww] = chebyshevPoints(nColPts,[-1,1]);
cheb.xx = xx;
cheb.ww = ww;
cheb.nSegment = nSegment;
cheb.nColPts = nColPts;

% Interpolate the guess at the chebyshev-points for transcription:
guess.tSpan = G.time([1,end]);
guess.time = getMultiChebTime(cheb,guess.tSpan);
guess.state = interp1(G.time', G.state', guess.time')';
guess.control = interp1(G.time', G.control', guess.time')';

[zGuess, pack] = packDecVar(guess.time, guess.state, guess.control);

% Unpack all bounds:
nGrid = nSegment*nColPts;
tLow = getMultiChebTime(cheb,[B.initialTime.low, B.finalTime.low]);
xLow = [B.initialState.low, B.state.low*ones(1,nGrid-2), B.finalState.low];
uLow = B.control.low*ones(1,nGrid);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = getMultiChebTime(cheb,[B.initialTime.upp, B.finalTime.upp]);
xUpp = [B.initialState.upp, B.state.upp*ones(1,nGrid-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nGrid);
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
nlpTime = toc;

soln = formatTrajectory(zSoln, pack, cheb);

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

[t, w] = getMultiChebTime(cheb,[z(1),z(2)]);
x = reshape(z((2+1):(2+nx)),nState,nTime);
u = reshape(z((2+nx+1):(2+nx+nu)),nControl,nTime);


end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [time, weights] = getMultiChebTime(cheb,tSpan)
%
% This function computes the time grid for a trajectory that is made up of
% a series of chebyshev orthogonal polynomials.
%
% INPUTS:
%   cheb = struct of information about the chebyshev basis functions
%       .xx = chebyshev points, on the domain [-1,1]
%       .ww = chebyshev weights, on the domain [-1,1]
%       .nSegment = number of trajectory segments
%   tSpan = [tInitial, tFinal]
%
% OUTPUTS:
%   time = [timeSegment_1, timeSegment_2, ... timeSegment_nSegment];
%
% NOTES:
%   This function will return duplicate times at the boundary to each
%   segment, so that the dynamics of each segment can easily be solved
%   independantly. These redundant points must be removed before returning
%   the output to the user.
%
%   For example, time should like something like:
%       time = [0, 1, 2, 3,  3, 4, 5, 6,  6, 7, 8, 9];
%

d = [0, (tSpan(2)-tSpan(1))/cheb.nSegment];   %Domain for the scaled points
[x, w] = chebyshevScalePoints(cheb.xx,cheb.ww,d);  %Scaled points

offset = linspace(tSpan(1),tSpan(2),cheb.nSegment+1);   %Starting time for each segment
time = x'*ones(1,cheb.nSegment) + ones(cheb.nColPts,1)*offset(1:(end-1));
nGrid = numel(time);
time = reshape(time,1,nGrid);

weights = reshape(w'*ones(1,cheb.nSegment),1,nGrid);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


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

nSegment = cheb.nSegment;

%%%% Enforce the dynamics:

% Analytic differentiation of the trajectory at chebyshev points:
domain = [0, (t(end)-t(1))/nSegment];   %Domain for the scaled points
xTemplate = chebyshevScalePoints(cheb.xx,cheb.ww,domain);  %Scaled points
D = chebyshevDifferentiationMatrix(xTemplate);
dxFun = zeros(size(x));
idx = 1:cheb.nColPts;
for i=1:nSegment  %Loop over each segment of the trajectory
    dxFun(:,idx) = (D*x(:,idx)')';
    idx = idx + cheb.nColPts;
end

% Derivative, according to the dynamics function:
dxDyn = dynFun(t,x,u);

% Add a constraint that both versions of the derivative must match.
% This ensures that the dynamics inside of each segment are correct.
dxError = dxFun - dxDyn;

% Add an additional defect that makes the state at the end of one
% segment match the state at the beginning of the next segment.
idxLow = cheb.nColPts*(1:(nSegment-1));
idxUpp = idxLow + 1;
stitchState = x(:,idxLow)-x(:,idxUpp);
stitchControl = u(:,idxLow)-u(:,idxUpp);  %Also need continuous control

defects = [...
    reshape(dxError, numel(dxError),1);
    reshape(stitchState, numel(stitchState),1);
    reshape(stitchControl, numel(stitchControl),1)];

%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u, defects, pathCst, bndCst);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function soln = formatTrajectory(zSoln, pack, cheb)
%
% This function formats the result of the trajectory optimization so that
% it is easy to use for plotting and analysis by the user.
%

[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack,cheb);

nSegment = cheb.nSegment;
nColPts = cheb.nColPts;

%%%% Need to remove the duplicate data points between segments:
idxLow = nColPts*(1:(nSegment-1));
idxUpp = idxLow + 1;
tSoln(idxUpp) = [];
xSoln(:,idxLow) = 0.5*(xSoln(:,idxLow)+xSoln(:,idxUpp));
xSoln(:,idxUpp) = [];
uSoln(:,idxLow) = 0.5*(uSoln(:,idxLow)+uSoln(:,idxUpp));
uSoln(:,idxUpp) = [];

%%%% Store the results:
soln.grid.time = tSoln;
soln.grid.state = xSoln;
soln.grid.control = uSoln;


%%%% Set up interpolating of the chebyshev trajectory:
idxKnot = (nColPts-1)*(1:(nSegment-1)) + 1;
soln.interp.state = @(t)( chebyshevMultiInterpolate(xSoln,tSoln,idxKnot,t) );
soln.interp.control = @(t)( chebyshevMultiInterpolate(uSoln,tSoln,idxKnot,t) );

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


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


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


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


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function D = chebyshevDifferentiationMatrix(x)
%
% Computes the chebyshev differentiation matrix
%
% NOTES:
%
% Example usage:   Df = (D*f')';
%       where  f = [nState x nPoints] values at each chebyshev node
%

n = length(x);

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

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

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

D = chebyshevDifferentiationMatrix(x);

%Apply the differentiation matrix
Df = (D*f')';

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function y = chebyshevMultiInterpolate(yData,tData,idxKnot,t)
%
% This function is a wrapper for chebyshevInterpolate that handles the
% piece-wise chebyshev polynomial trajectory
%

% All points are considered valid, so extend edge bins:
Tbins = [-inf, tData(idxKnot), inf];

%Figure out which bins each query is in:
[~, idx] = histc(t,Tbins);

% Loop over each segment of the trajectory:
ny = size(yData,1);
nt = length(t);
gridIdx = [1, idxKnot, length(tData)];
nSegment = length(gridIdx)-1;
y = zeros(ny,nt);
for i=1:nSegment
    if sum(idx==i)>0 % Then there are points to evaluate here!
        y(:,(idx==i)) = chebyshevInterpolate(...
            yData(:,gridIdx(i):gridIdx(i+1)),...
            t(idx==i),...
            tData(gridIdx([i,i+1]))  );
    end
end

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


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


