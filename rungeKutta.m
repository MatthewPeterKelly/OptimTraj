function soln = rungeKutta(problem,defaultOptions)
% soln = rungeKutta(problem)
%
% This function transcribes a trajectory optimization problem using the
% multiple shooting, with 4th-order Runge Kutta integration
%
% See Bett's book for details on the method
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'rungeKutta'
%   problem.options.rungeKutta = struct with method parameters:
%       .nSegment = number of trajectory segments
%       .nSubStep = number of sub-steps to use in each segment
%


%%%%TODO%%%%
% 
% 1) create locally global variables to prevent double-calling the dynamics function
% 
% 2) create a nice wrapper for RK4 that computes all intermediate states, for constraints and cost function
%
% 3) Use RK4 to compute the cost function as well as the defects. It makes sense, just think a bit harder about it.
%
%%%%%%%%%%%%




%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

% Figure out grid size:
nSegment = Opt.rungeKutta.nSegment;
nSubStep = Opt.rungeKutta.nSubStep;
nControl = 2*nSegment*nSubStep + 1;
nState = nSegment + 1;

% Print out some solver info if desired:
if Opt.verbose > 0
    fprintf('  -> Transcription via 4th-order Runge-Kutta method \n');
    fprintf('        nSegments = %d \n', nSegment);
    fprintf('        nSubSteps = %d \n', nSubStep);
end

% Interpolate the guess at the grid-points for transcription:
guess.tSpan = G.time([1,end]);
guess.state = interp1(G.time', G.state', guess.time')';
guess.control = interp1(G.time', G.control', guess.time')';

[zGuess, pack] = packDecVar(guess.time, guess.state, guess.control);

% Unpack all bounds:
tLow = linspace(B.initialTime.low, B.finalTime.low, nGrid);
xLow = [B.initialState.low, B.state.low*ones(1,nGrid-2), B.finalState.low];
uLow = B.control.low*ones(1,nGrid);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = linspace(B.initialTime.upp, B.finalTime.upp, nGrid);
xUpp = [B.initialState.upp, B.state.upp*ones(1,nGrid-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nGrid);
zUpp = packDecVar(tUpp,xUpp,uUpp);

%%%% Set up problem for fmincon:

P.objective = @(z)( ...
    myObjective(z, pack, F.pathObj, F.bndObj) );

P.nonlcon = @(z)( ...
    myConstraint(z, pack, F.dynamics, F.pathCst, F.bndCst) );

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
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack);
nlpTime = toc;

%%%% Store the results:

soln.grid.time = tSoln;
soln.grid.state = xSoln;
soln.grid.control = uSoln;

% Use quadratic interpolation for each trajectory segment
soln.interp.state = @(t)( pwPoly2(tSoln,xSoln,t) );
soln.interp.control = @(t)( pwPoly2(tSoln,uSoln,t) );

soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                   SUB FUNCTIONS                                   %%%%
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

function [t,x,u] = unPackDecVar(z,pack)
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
%

nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;
nx = nState*nTime;
nu = nControl*nTime;

t = linspace(z(1),z(2),nTime);

x = reshape(z((2+1):(2+nx)),nState,nTime);
u = reshape(z((2+nx+1):(2+nx+nu)),nControl,nTime);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function cost = myObjective(z,pack,pathObj,bndObj)
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

[t,x,u] = unPackDecVar(z,pack);

% Compute the cost integral along trajectory
if isempty(pathObj)
    integralCost = 0;
else
    dt = 2*(t(end)-t(1))/(pack.nTime-1);
    integrand = pathObj(t,x,u);  %Calculate the integrand of the cost function
    
    %Simpson quadrature for integration of the cost function:
    weights = 2*ones(pack.nTime,1);
    weights(2:2:end) = 4;
    weights([1,end]) = 1;
    weights = weights/6;
    
    integralCost = dt*integrand*weights;  %Trapazoidal integration
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

function [c, ceq] = myConstraint(z,pack,dynFun, pathCst, bndCst)
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

[t,x,u] = unPackDecVar(z,pack);


%%%% Compute dynamics at each grid point:
dt = 2*(t(end)-t(1))/(pack.nTime-1);
dx = dynFun(t,x,u);

iLow = 1:2:(pack.nTime-1);
iMid = iLow + 1;
iUpp = iMid + 1;

xLow = x(:,iLow);
xMid = x(:,iMid);
xUpp = x(:,iUpp);

fLow = dx(:,iLow);
fMid = dx(:,iMid);
fUpp = dx(:,iUpp);


%%%% Mid-point constraint (Hermite)
defectMidpoint = xMid - (xUpp+xLow)/2 - dt*(fLow-fUpp)/8;

%%%% Interval constraint (Simpson)
defectInterval = xUpp - xLow - dt*(fUpp + 4*fMid + fLow)/6;


%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,...
    [defectMidpoint,defectInterval],...
    pathCst, bndCst);

end