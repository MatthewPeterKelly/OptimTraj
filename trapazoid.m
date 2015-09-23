function soln = trapazoid(problem)
% soln = trapazoid(problem)
%
% This function transcribes a trajectory optimization problem using the
% trapazoid method for enforcing the dynamics. It can be found in chapter
% four of Bett's book:
%
%   John T. Betts, 2001
%   Practical Methods for Optimal Control Using Nonlinear Programming
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'trapazoid'
%   problem.options.trapazoid = struct with method parameters:
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
%                 dObjGrad = [1, 1+nx+nu, nTime]
%
%         [obj, objGrad] = bndObj(t0,x0,tF,xF)
%                 obj = scalar = objective function for boundry points
%                 objGrad = [1, 1+nx+1+nx]
%
%         [c, ceq, cGrad, ceqGrad] = pathCst(t,x,u)
%                 c = column vector of inequality constraints  ( c <= 0 )
%                 ceq = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+nu, nTime];
%                 ceqGrad = [nCstEq, 1+nx+nu, nTime];
%
%         [c, ceq, cGrad, ceqGrad] = bndCst(t0,x0,tF,xF)
%                 c = column vector of inequality constraints  ( c <= 0 )
%                 ceq = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+1+nx];
%                 ceqGrad = [nCstEq, 1+nx+1+nx];
%
% The gradients should be taken with respect to the inputs
%

%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

nGrid = Opt.trapazoid.nGrid;  %Number of grid points for transcription

flagGradObj = strcmp(Opt.nlpOpt.GradObj,'on');
flagGradCst = strcmp(Opt.nlpOpt.GradConstr,'on');

% Print out some solver info if desired:
if Opt.verbose > 0
    fprintf('  -> Transcription via trapazoid method, nGrid = %d\n ',Opt.trapazoid.nGrid);
    if flagGradObj
        fprintf('      - using analytic gradients of objective function\n');
    end
    if flagGradCst
        fprintf('      - using analytic gradients of constraint function\n');
    end
    fprintf('\n');
end

% Interpolate the guess at the grid-points for transcription:
guess.tSpan = G.time([1,end]);
guess.time = linspace(guess.tSpan(1), guess.tSpan(2), nGrid);
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
if flagGradObj
    P.objective = @(z)( ...
        myObjGrad(z, pack, F.pathObj, F.bndObj) );   %Analytic gradients
else
    P.objective = @(z)( ...
        myObjective(z, pack, F.pathObj, F.bndObj) );   %Numerical gradients
end
if flagGradCst
    P.nonlcon = @(z)( ...
        myCstGrad(z, pack, F.dynamics, F.pathCst, F.bndCst) ); %Analytic gradients
else
    P.nonlcon = @(z)( ...
        myConstraint(z, pack, F.dynamics, F.pathCst, F.bndCst) ); %Numerical gradients
end
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

soln.interp.state = @(t)( interp1(tSoln',xSoln',t','linear',nan)' );
soln.interp.control = @(t)( interp1(tSoln',uSoln',t','linear',nan)' );

soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
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
    dt = (t(end)-t(1))/(pack.nTime-1);
    integrand = pathObj(t,x,u);  %Calculate the integrand of the cost function
    weights = ones(length(t),1); weights([1,end]) = 0.5;
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


%%%% Compute defects along the trajectory:

dt = (t(end)-t(1))/(pack.nTime-1);
dx = dynFun(t,x,u);

xLow = x(:,1:end-1);
xUpp = x(:,2:end);

dxLow = dx(:,1:end-1);
dxUpp = dx(:,2:end);

% This is the key line:  (Trapazoid Rule)
defects = xUpp-xLow - 0.5*dt*(dxLow+dxUpp);


%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,defects, pathCst, bndCst);

end




%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%            Additional Sub-Functions for Gradients                 %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function [dt, dtGrad] = getTimeStepGrad(t,tIdx,nz)
%
% Computes the time step and its gradient
%
% dt = [1,1]
% dtGrad = [1,nz]
%

nTime = length(t);

dt = (t(end)-t(1))/(nTime-1);
dtGrad = zeros(1,nz);

dtGrad(1,tIdx(1)) = -1/(nTime-1);
dtGrad(1,tIdx(end)) = 1/(nTime-1);

end



function [tGrad, xGrad, uGrad] = getDecVarGrad(t,x,u,tIdx,xIdx,uIdx,nz)
%
% Computes the gradients of t,x,u with respect to the original decision
% variable vector.
%

nt = size(t,2);
nx = size(x,1);
nu = size(u,1);

tGrad = zeros(1,nt,nz);
xGrad = zeros(nx,nt,nz);
uGrad = zeros(nu,nt,nz);

%%%% Compute gradients of time:
% alpha = (0..N-1)/(N-1)
% t = alpha*tUpp + (1-alpha)*tLow
alpha = (0:(nt-1))/(nt-1);
tGrad(1,:,tIdx(1)) = 1-alpha;
tGrad(1,:,tIdx(2)) = alpha;

%%%% Compute gradients of state and control:
for i=1:nt
    for j=1:nx
        xGrad(j,i,xIdx(j,i)) = 1;
    end
    for j=1:nu
        uGrad(j,i,uIdx(j,i)) = 1;
    end
end

end


%%%%%%%%%%%%%%%%%


function [dCost, dCostGrad] = getObjGrad(t,x,u,tGrad, xGrad, uGrad,pathObj)


[dCost, integrandGrad] = pathObj(t,x,u);  %Calculate the integrand of the cost function

dCost = dCost';

stateGrad = [tGrad;xGrad;uGrad];

[~,nt,nz] = size(stateGrad);
dCostGrad = zeros(nt,nz);

for i=1:nt
    for j=1:nz
        dCostGrad(i,j) = sum(integrandGrad(:,i).*stateGrad(:,i,j));
    end
end

end


%%%%%%%%%%%%%%%%%


function [cost, costGrad] = myObjGrad(z,pack,pathObj,bndObj)
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

% Dummy vector for tracking indicies:
nz = length(z);
index = 1:nz;
[tIdx,xIdx,uIdx] = unPackDecVar(index,pack);  tIdx = tIdx([1,end]);

%Unpack the decision variables:
[t,x,u] = unPackDecVar(z,pack);

    % Time step for integration:
    [dt, dtGrad] = getTimeStepGrad(t,tIdx,nz);
nx = size(x,1);
nu = size(u,1);

% Compute the cost integral along the trajectory
if isempty(pathObj)
    integralCost = 0;
    integralCostGrad = zeros(nx+nu,1);
else
    
    
    %Compute gradients of the decision variables:
    [tGrad, xGrad, uGrad] = getDecVarGrad(t,x,u,tIdx,xIdx,uIdx,nz);
    
    % integration weights:
    weights = ones(pack.nTime,1); weights([1,end]) = 0.5;
    
    % Objective function integrand and gradients:
    [dCost, dCostGrad] = getObjGrad(t,x,u,tGrad, xGrad, uGrad,pathObj);
    
    % integral objective function
    integralCost = dt*weights'*dCost;
    
    % Gradient of integral objective function
    integralCostGrad = ...
        dtGrad*(weights'*dCost) + ...
        dt*weights'*dCostGrad;
end

% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
    bndCostGrad = zeros(1,nz);
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    [bndCost, bndCostGrad] = bndObj(t0,x0,tF,xF);
    error('Not yet implemented!');
end

cost = bndCost + integralCost;

% Gradients...
costGrad = bndCostGrad + integralCostGrad;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq, cGrad, ceqGrad] = myCstGrad(z,pack,dynFun, pathCst, bndCst)
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


%%%% Compute defects along the trajectory:

dt = (t(end)-t(1))/(pack.nTime-1);
[dx, dxGrad] = dynFun(t,x,u);

xLow = x(:,1:end-1);
xUpp = x(:,2:end);

dxLow = dx(:,1:end-1);
dxUpp = dx(:,2:end);

% This is the key line:  (Trapazoid Rule)
defects = xUpp-xLow - 0.5*dt*(dxLow+dxUpp);


%%%% Call user-defined constraints and pack up:

ceq_dyn = reshape(defects,numel(defects),1);

%%%% Compute the user-defined constraints:
if isempty(pathCst)
    c_path = [];
    ceq_path = [];
else
    [c_path, ceq_path, c_pathGrad, ceq_pathGrad] = pathCst(t,x,u);
end
if isempty(bndCst)
    c_bnd = [];
    ceq_bnd = [];
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    [c_bnd, ceq_bnd, c_bndGrad, ceq_bndGrad] = bndCst(t0,x0,tF,xF);
end

%%%% Pack everything up:
c = [c_path;c_bnd];
ceq = [ceq_dyn; ceq_path; ceq_bnd];

%%%% TODO:
cGrad = [];
ceqGrad = [];

end

