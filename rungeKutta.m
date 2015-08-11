function soln = rungeKutta(problem)
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

%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

% Figure out grid size:
nSegment = Opt.rungeKutta.nSegment;
nSubStep = Opt.rungeKutta.nSubStep;
nGridControl = 2*nSegment*nSubStep + 1;
nGridState = nSegment + 1;

% Print out some solver info if desired:
if Opt.verbose > 0
    fprintf('  -> Transcription via 4th-order Runge-Kutta method \n');
    fprintf('        nSegments = %d \n', nSegment);
    fprintf('        nSubSteps = %d \n', nSubStep);
end

% Interpolate the guess at the transcription grid points for initial guess:
guess.tSpan = G.time([1,end]);
guess.tState = linspace(guess.tSpan(1), guess.tSpan(2), nGridState);
guess.tControl = linspace(guess.tSpan(1), guess.tSpan(2), nGridControl);
guess.state = interp1(G.time', G.state', guess.tState')';
guess.control = interp1(G.time', G.control', guess.tControl')';
[zGuess, pack] = packDecVar(guess.tSpan, guess.state, guess.control);

% Unpack all bounds:
tLow = [B.initialTime.low, B.finalTime.low];
xLow = [B.initialState.low, B.state.low*ones(1,nGridState-2), B.finalState.low];
uLow = B.control.low*ones(1,nGridControl);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = [B.initialTime.upp, B.finalTime.upp];
xUpp = [B.initialState.upp, B.state.upp*ones(1,nGridState-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nGridControl);
zUpp = packDecVar(tUpp,xUpp,uUpp);

%%%% Set up problem for fmincon:

P.objective = @(z)( ...
    myObjective(z, pack, F.dynamics, F.pathObj, F.bndObj, F.pathCst, F.bndCst) );

P.nonlcon = @(z)( myConstraint(z) );

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
[tSpan,~,uSoln] = unPackDecVar(zSoln,pack);
nlpTime = toc;

%%%% Store the results:
[~,tGrid,xGrid,uGrid] = P.objective(zSoln);
soln.grid.time = tGrid;
soln.grid.state = xGrid;
soln.grid.control = uGrid;

% Quadratic interpolation over each sub-step for the control:
tSoln = linspace(tSpan(1),tSpan(2),nGridControl);
soln.interp.control = @(t)( pwPoly2(tSoln, uSoln, t) );

% Cubic spline representation of the state over each substep:
dxGrid = F.dynamics(tGrid,xGrid,uGrid);
xSpline = pwch(tGrid, xGrid, dxGrid);
soln.interp.state = @(t)( ppval(xSpline,t) );

% General information about the optimization run
soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                   SUB FUNCTIONS                                   %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [decVars,pack] = packDecVar(tSpan,state,control)
%
% This function collapses the time (t), state (x)
% and control (u) matricies into a single vector
%
% INPUTS:
%   tSpan = [1, 2] = time bounds
%   state = [nState, nGridState] = state vector at each grid point
%   control = [nControl, nGridControl] = control vector at each grid point
%
% OUTPUTS:
%   decVars = column vector of 2 + nState*nGridState + nControl*nGridControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nState
%       .nGridState
%       .nControl
%       .nGridControl
%
% NOTES:
% nGridControl = 2*nSegment*nSubStep + 1;
% nGridState = nSegment + 1;
%

[nState, nGridState] = size(state);
[nControl, nGridControl] = size(control);

xCol = reshape(state, nState*nGridState, 1);
uCol = reshape(control, nControl*nGridControl, 1);

decVars = [tSpan(:);xCol;uCol];

pack.nState = nState;
pack.nGridState = nGridState;
pack.nControl = nControl;
pack.nGridControl = nGridControl;
pack.nSegment = nGridState - 1;
pack.nSubStep = (nGridControl-1)/(2*pack.nSegment);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [tSpan, state, control] = unPackDecVar(decVars,pack)
%
% This function unpacks the decision variables for
% trajectory optimization into the time (t),
% state (x), and control (u) matricies
%
% INPUTS:
%   decVars = column vector of 2 + nState*nGridState + nControl*nGridControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nState
%       .nGridState
%       .nControl
%       .nGridControl
%
% OUTPUTS:
%   tSpan = [1, 2] = time bounds
%   state = [nState, nGridState] = state vector at each grid point
%   control = [nControl, nGridControl] = control vector at each grid point
%

nx = pack.nState*pack.nGridState;
nu = pack.nControl*pack.nGridControl;

tSpan = [decVars(1),decVars(2)];

state = reshape(decVars((2+1):(2+nx)), pack.nState, pack.nGridState);
control = reshape(decVars((2+nx+1):(2+nx+nu)), pack.nControl, pack.nGridControl);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [cost, t,x,u] = myObjective(decVars, pack, ...
    dynamics, pathObj, bndObj, pathCst, bndCst)
%
% This function unpacks the decision variables, sends them to the
% user-defined objective functions, and then returns the final cost
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   all = user-defined integral objective function
%   endObj = user-defined end-point objective function
%
% OUTPUTS:
%   cost = scale cost for this set of decision variables
%
% NOTES:
%   All calculations for both the objective function and the constraints
%   are done here, to prevent double-calling the dynamics function. FMINCON
%   always calls the obejctive function and then the constraint function
%   with the same state vector, so both are calculated here, and the
%   constraints are stored as a global variable to be read later by the
%   constraint function.
%

global RUNGE_KUTTA_CONSTRAINT_CEQ
global RUNGE_KUTTA_CONSTRAINT_C
global RUNGE_KUTTA_CHECK_DEC_VARS
RUNGE_KUTTA_CHECK_DEC_VARS = decVars;

[tSpan, state, control] = unPackDecVar(decVars,pack);

nState = pack.nState;
nControl = pack.nControl;
nSegment = pack.nSegment;
nSubStep = pack.nSubStep;

% time, state, and control at the ends of each substep
nTime = 1+nSegment*nSubStep;
t = linspace(tSpan(1), tSpan(2), nTime);
x = zeros(nState, nTime);
u = control(nControl,1:2:end);  %Omit the control at the midpoint
c = zeros(1, nTime-1);  %Integral cost for each segment
dt = (t(end)-t(1))/(nTime-1);

idx = 1:nSubStep:(nTime-1);   %Indicies for the start of each segment
x(:,[idx,end]) = state;   %Fill in the states that we already know
for iSubStep = 1:nSubStep
% March forward Runge-Kutta step   

t0 = t(idx);
x0 = x(:,idx);
uLow = u(:,idx);
uMid = control(:,2*idx);
uUpp = u(:,idx+1);
k0 = combinedDynamics(t0,        x0,                         uLow, dynamics,pathObj);
k1 = combinedDynamics(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid, dynamics,pathObj);
k2 = combinedDynamics(t0+0.5*dt, x0 + 0.5*dt*k1(1:nState,:), uMid, dynamics,pathObj);
k3 = combinedDynamics(t0+dt,     x0 +     dt*k2(1:nState,:), uUpp, dynamics,pathObj);
z = (dt/6)*(k0 + 2*k1 + 2*k2 + k3);  %Change over the sub-step
xNext = x0 + z(1:nState,:);  %Next state
c(idx) = z(end,:);  %Integral of the cost function over this step

if iSubStep == nSubStep %We've reached the end of the interval
    % Compute the defect vector:
    defects = xNext - x(:,idx+1);
else
    % Store the state for next step in time
    idx = idx+1;
    x(:,idx) = xNext;
end

end

% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
else
    t0 = tSpan(1);
    tF = tSpan(2);
    x0 = state(:,1);
    xF = state(:,end);
    bndCost = bndObj(t0,x0,tF,xF);
end

integralCost = sum(c);  %Sum up the integral cost over each segment
cost = bndCost + integralCost;

%%%% Call user-defined constraints and pack up:
[RUNGE_KUTTA_CONSTRAINT_C, RUNGE_KUTTA_CONSTRAINT_CEQ] = ...
    collectConstraints(t,x,u,...
    defects,...
    pathCst, bndCst);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq] = myConstraint(decVars)
%
% This function unpacks the decision variables, computes the defects along
% the trajectory, and then evaluates the user-defined constraint functions.
%
global RUNGE_KUTTA_CONSTRAINT_CEQ
global RUNGE_KUTTA_CONSTRAINT_C
global RUNGE_KUTTA_CHECK_DEC_VARS

matchFail = decVars ~= RUNGE_KUTTA_CHECK_DEC_VARS;
if any(matchFail)
    disp('ERROR: Match failed!')
end

c = RUNGE_KUTTA_CONSTRAINT_C;
ceq = RUNGE_KUTTA_CONSTRAINT_CEQ;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function dz = combinedDynamics(t,x,u,dynamics,pathObj)
% dz = combinedDynamics(t,x,u,dynamics,pathObj)
%
% This function packages the dynamics and the cost function together so
% that they can be integrated at the same time. 
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%   dynamics(t,x,u) = dynamics function handle
%               dx = [nState, nTime] = dx/dt = derivative of state wrt time
%   pathObj(t,x,u) = integral cost function handle
%                 dObj = [1, nTime] = integrand from the cost function
%
% OUTPUTS:
%   dz = [dx; dObj] = combined dynamics of state and cost
%

dx = dynamics(t,x,u);
if isempty(pathObj)
    dc = zeros(size(t));
else
    dc = pathObj(t,x,u);
end

dz = [dx;dc];  %Combine and return

end
