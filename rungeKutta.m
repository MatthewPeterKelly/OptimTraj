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
%       .AdaptiveDerivativeCheck = 'off' by default. Set to 'on' to enable
%           numerical checks on the analytic gradients, computed using the
%           derivest package, rather than fmincon's internal checks.
%           Derivest is slower, but more accurate than fmincon. Derivest
%           can be downloaded from the Mathworks File Exchange, file id of
%           13490 - Adaptive Robust Numerical Differentation, John D-Errico
%
%
% NOTES:
%
%   Code for computing analyic gradients of the Runge Kutta method was
%   contributed by Will Wehner.
%
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
flagGradObj = strcmp(Opt.nlpOpt.GradObj,'on');
flagGradCst = strcmp(Opt.nlpOpt.GradConstr,'on');

% gradient info for use in calculating analytic gradients of objective and
% constraints
if flagGradCst || flagGradObj
    gradInfo = grad_computeInfo(pack);
    disp('  -> Using analytic gradients');
end

if flagGradObj
    P.objective = @(z)( ...
        myObjGrad(z, pack, F.dynamics, F.pathObj, F.bndObj, gradInfo) );   %Analytic gradients
else
    P.objective = @(z)( ...
        myObjective(z, pack, F.dynamics, F.pathObj, F.bndObj) );   %Numerical gradients
end

if flagGradCst
    P.nonlcon = @(z)( ...
        myCstGrad(z, pack, F.dynamics, F.pathObj, F.pathCst, F.bndCst, gradInfo) ); %Analytic gradients
else
    P.nonlcon = @(z)( ...
        myConstraint(z, pack, F.dynamics, F.pathObj, F.pathCst, F.bndCst) ); %Numerical gradients
end

% Check analytic gradients with DERIVEST package
if strcmp(Opt.rungeKutta.AdaptiveDerivativeCheck,'on')
    if exist('jacobianest','file')
        runGradientCheck(zGuess, pack,F.dynamics, F.pathObj, F.bndObj, F.pathCst, F.bndCst, gradInfo);
        Opt.nlpOpt.DerivativeCheck = [];  %Disable built-in derivative check
    else
        Opt.rungeKutta.AdaptiveDerivativeCheck = 'cannot find jacobianest.m';
        disp('Warning: the derivest package is not on search path.');
        disp(' --> Using fmincon''s built-in derivative checks.');
    end
end

% Build the standard fmincon problem struct
P.x0 = zGuess;
P.lb = zLow;
P.ub = zUpp;
P.Aineq = []; P.bineq = [];
P.Aeq = []; P.beq = [];
P.solver = 'fmincon';
P.options = Opt.nlpOpt;

%%%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal,exitFlag,output] = fmincon(P);
[tSpan,~,uSoln] = unPackDecVar(zSoln,pack);
nlpTime = toc;

%%%% Store the results:
[tGrid,xGrid,uGrid] = simulateSystem(zSoln, pack, F.dynamics, F.pathObj);
soln.grid.time = tGrid;
soln.grid.state = xGrid;
soln.grid.control = uGrid;

% Quadratic interpolation over each sub-step for the control:
tSoln = linspace(tSpan(1),tSpan(2),nGridControl);
soln.interp.control = @(t)( interp1(tSoln', uSoln', t','pchip')' );

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

function cost = myObjective(decVars, pack,dynamics, pathObj, bndObj)
%
% This function unpacks the decision variables, sends them to the
% user-defined objective functions, and then returns the final cost
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   bndObj = user-defined boundary objective function
%
% OUTPUTS:
%   cost = scalar cost for this set of decision variables
%
%

% All of the real work happens inside this function:
[t,x,~,~,pathCost] = simulateSystem(decVars, pack, dynamics, pathObj);

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

cost = bndCost + pathCost;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq] = myConstraint(decVars, pack, dynamics, pathObj, pathCst, bndCst)
%
% This function unpacks the decision variables, computes the defects along
% the trajectory, and then evaluates the user-defined constraint functions.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   pathCst = user-defined path-constraint function
%   bndCst = user-defined boundary constraint function
%
% OUTPUTS:
%   c = non-linear inequality constraint
%   ceq = non-linear equatlity cosntraint
%
% NOTE:
%   - path constraints are  satisfied at the start and end of each sub-step
%


[t,x,u,defects] = simulateSystem(decVars, pack, dynamics, pathObj);

%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,...
    defects,...
    pathCst, bndCst);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [t,x,u,defects,pathCost] = simulateSystem(decVars, pack, dynFun, pathObj)
%
% This function does the real work of the transcription method. It
% simulates the system forward in time across each segment of the
% trajectory, computes the integral of the cost function, and then matches
% up the defects between the end of each segment and the start of the next.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%
% OUTPUTS:
%   t = [1 x nGrid] = time vector for the edges of the sub-step grid
%   x = [nState x nGrid] = state vector
%   u = [nControl x nGrid] = control vector
%   defects = [nState x nSegment] = defect matrix
%   pathCost = scalar cost for the path integral
%
% NOTES:
%   - nGrid = nSegment*nSubStep+1
%   - This function is usually called twice for each combination of
%   decision variables: once by the objective function and once by the
%   constraint function. To keep the code fast I cache the old values and
%   only recompute when the inputs change.
%


%%%% CODE OPTIMIZATION %%%%
%
% Prevents the same exact code from being called twice by caching the
% solution and reusing it when appropriate.
%
global RUNGE_KUTTA_t RUNGE_KUTTA_x RUNGE_KUTTA_u
global RUNGE_KUTTA_defects RUNGE_KUTTA_pathCost
global RUNGE_KUTTA_decVars
%
usePreviousValues = false;
if ~isempty(RUNGE_KUTTA_decVars)
    if length(RUNGE_KUTTA_decVars) == length(decVars)
        if ~any(RUNGE_KUTTA_decVars ~= decVars)
            usePreviousValues = true;
        end
    end
end
%
if usePreviousValues
    t = RUNGE_KUTTA_t;
    x = RUNGE_KUTTA_x;
    u = RUNGE_KUTTA_u;
    defects = RUNGE_KUTTA_defects;
    pathCost = RUNGE_KUTTA_pathCost;
else
    %
    %
    %%%% END CODE OPTIMIZATION %%%%
    
    
    [tSpan, state, control] = unPackDecVar(decVars,pack);
    
    nState = pack.nState;
    nSegment = pack.nSegment;
    nSubStep = pack.nSubStep;
    
    % NOTES:
    %   The following bit of code is a bit confusing, mostly due to the
    %   need for vectorization to make things run at a reasonable speed in
    %   Matlab. Part of the confusion comes because the decision variables
    %   include the state at the beginning of each segment, but the control
    %   at the beginning and middle of each substep - thus there are more
    %   control grid-points than state grid points. The calculations are
    %   vectorized over segments, but not sub-steps, since the result of
    %   one sub-step is required for the next.
    
    % time, state, and control at the ends of each substep
    nTime = 1+nSegment*nSubStep;
    t = linspace(tSpan(1), tSpan(2), nTime);
    x = zeros(nState, nTime);
    u = control(:,1:2:end); % Control a the endpoints of each segment
    uMid = control(:,2:2:end);  %Control at the mid-points of each segment
    c = zeros(1, nTime-1);  %Integral cost for each segment
    dt = (t(end)-t(1))/(nTime-1);
    
    idx = 1:nSubStep:(nTime-1);   %Indicies for the start of each segment
    x(:,[idx,end]) = state;   %Fill in the states that we already know
    
    for iSubStep = 1:nSubStep
        % March forward Runge-Kutta step
        
        t0 = t(idx);
        x0 = x(:,idx);
        
        k0 = combinedDynamics(t0,        x0,                         u(:,idx), dynFun,pathObj);
        k1 = combinedDynamics(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid(:,idx), dynFun,pathObj);
        k2 = combinedDynamics(t0+0.5*dt, x0 + 0.5*dt*k1(1:nState,:), uMid(:,idx), dynFun,pathObj);
        k3 = combinedDynamics(t0+dt,     x0 +     dt*k2(1:nState,:), u(:,idx+1), dynFun,pathObj);
        z = (dt/6)*(k0 + 2*k1 + 2*k2 + k3);  %Change over the sub-step
        
        xNext = x0 + z(1:nState,:);  %Next state
        c(idx) = z(end,:);  %Integral of the cost function over this step
        
        if iSubStep == nSubStep %We've reached the end of the interval
            % Compute the defect vector:
            defects = xNext - x(:,idx+1);
        else
            % Store the state for next step in time
            idx = idx+1;   %  <-- This is important!!
            x(:,idx) = xNext;
        end
        
    end
    
    pathCost = sum(c);  %Sum up the integral cost over each segment
    
    %%%% Cache results to use on the next call to this function.
    RUNGE_KUTTA_t = t;
    RUNGE_KUTTA_x = x;
    RUNGE_KUTTA_u = u;
    RUNGE_KUTTA_defects = defects;
    RUNGE_KUTTA_pathCost = pathCost;
    
end

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function dz = combinedDynamics(t,x,u,dynFun,pathObj)
% dz = combinedDynamics(t,x,u,dynFun,pathObj)
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


dx = dynFun(t,x,u);
if isempty(pathObj)
    dc = zeros(size(t));
else
    dc = pathObj(t,x,u);
end

dz = [dx;dc];  %Combine and return


end













%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                 Analytic Gradient Stuff                           %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

















function gradInfo = grad_computeInfo(pack)
%
% This function computes the matrix dimensions and indicies that are used
% to map the gradients from the user functions to the gradients needed by
% fmincon. The key difference is that the gradients in the user functions
% are with respect to their input (t,x,u) or (t0,x0,tF,xF), while the
% gradients for fmincon are with respect to all decision variables.
%
% INPUTS:
%   nDeVar = number of decision variables
%   pack = details about packing and unpacking the decision variables
%       .nTime
%       .nState
%       .nControl
%
% OUTPUTS:
%   gradInfo = details about how to transform gradients
%


%nTime = pack.nTime;
nState = pack.nState;
nGridState = pack.nGridState;
nControl = pack.nControl;
nGridControl = pack.nGridControl;

nDecVar = 2 + nState*nGridState + nControl*nGridControl;

zIdx = 1:nDecVar;
gradInfo.nDecVar = nDecVar;
[tIdx, xIdx, uIdx] = unPackDecVar(zIdx,pack);
gradInfo.tIdx = tIdx([1,end]);
%gradInfo.xuIdx = [xIdx;uIdx];
gradInfo.xIdx = xIdx;
gradInfo.uIdx = uIdx;

nSegment = pack.nSegment;
nSubStep = pack.nSubStep;

% indices of decVars associated with u
indu = 1:2:(1+2*nSegment*nSubStep);
gradInfo.indu = uIdx(:,indu);
% indices of decVars associated with uMid
indumid = 2:2:(1+2*nSegment*nSubStep);
gradInfo.indumid = uIdx(:,indumid);

%%%% Compute gradients of time:
%%%% alpha = (0..N-1)/(N-1)
%%%% t = alpha*tUpp + (1-alpha)*tLow
% alpha = (0:(nTime-1))/(nTime-1);
% gradInfo.alpha = [1-alpha; alpha];
%
% if (gradInfo.tIdx(1)~=1 || gradInfo.tIdx(end)~=2)
%     error('The first two decision variables must be the initial and final time')
% end
% gradInfo.dtGrad = [-1; 1]/(nTime-1);

%%%% Compute gradients of state (What is this for)?
% gradInfo.xGrad = zeros(nState,nTime,nDecVar);
% for iTime=1:nTime
%     for iState=1:nState
%         gradInfo.xGrad(iState,iTime,xIdx(iState,iTime)) = 1;
%     end
% end

%%%% For unpacking the boundary constraints and objective:
gradInfo.bndIdxMap = [tIdx(1); xIdx(:,1); tIdx(end); xIdx(:,end)];


end



%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%




function [fail] = runGradientCheck(z_test, pack,dynamics, pathObj, bndObj, pathCst, bndCst, gradInfo)
%
% This function tests the analytic gradients of the objective and
% nonlinear constraints with the DERIVEST package. The finite difference
% calculations in matlab's optimization package were not sufficiently
% accurate.
%
GradientCheckTol = 1e-6;  %Analytic gradients must match numerical within this bound

fail = 0;

fprintf('\n%s\n','____________________________________________________________')
fprintf('%s\n','  DerivativeCheck Information with DERIVEST Package ')

% analytic gradient
[~, dcost] = myObjGrad(z_test, pack, dynamics, pathObj, bndObj, gradInfo);

% check gradient with derivest package
deriv = gradest(@(z) myObjGrad(z, pack, dynamics, pathObj, bndObj, gradInfo),z_test);

% print largest difference in numerical and analytic gradients
fprintf('\n%s\n','Objective function derivatives:')
fprintf('%s\n','Maximum relative difference between user-supplied')
fprintf('%s %1.5e \n','and finite-difference derivatives = ',max(abs(dcost-deriv')))
if any(abs(dcost-deriv') > GradientCheckTol)
    error('Objective gradient did not pass')
end

% analytic nonlinear constraints
[c, ceq,dc, dceq] = myCstGrad(z_test, pack, dynamics, pathObj, pathCst, bndCst, gradInfo);

% check nonlinear inequality constraints with 'jacobianest'
if ~isempty(c)
    jac = jacobianest(@(z) myConstraint(z, pack, dynamics, pathObj, pathCst, bndCst),z_test);
    
    % print largest difference in numerical and analytic gradients
    fprintf('\n%s\n','Nonlinear inequality constraint function derivatives:')
    fprintf('%s\n','Maximum relative difference between user-supplied')
    fprintf('%s %1.5e \n','and finite-difference derivatives = ',max(max(abs(dc-jac'))))
    if any(any(abs(dc - jac') > GradientCheckTol))
        error('Nonlinear inequality constraint did not pass')
    end
end

% check nonlinear equality constraints with 'jacobianest'
if ~isempty(ceq)
    jac = jacobianest(@(z) myCstGradCheckEq(z, pack, dynamics, pathObj, pathCst, bndCst),z_test);
    
    % print largest difference in numerical and analytic gradients
    fprintf('\n%s\n','Nonlinear equality constraint function derivatives:')
    fprintf('%s\n','Maximum relative difference between user-supplied')
    fprintf('%s %1.5e \n','and finite-difference derivatives = ',max(max(abs(dceq-jac'))))
    if any(any(abs(dceq - jac') > GradientCheckTol))
        error('Nonlinear equality constraint did not pass')
    end
end

fprintf('\n%s\n','DerivativeCheck successfully passed.')
fprintf('%s\n','____________________________________________________________')
end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function ceq = myCstGradCheckEq(decVars, pack, dynamics, pathObj, pathCst, bndCst)
% This function is necessary for runGradientCheck function
% return only equality constraint (ceq) for use with jacobest.m

[t,x,u,defects] = simulateSystem(decVars, pack, dynamics, pathObj);

%%%% Call user-defined constraints and pack up:
[~, ceq] = collectConstraints(t,x,u,...
    defects,...
    pathCst, bndCst);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%




function [cost, dcost] = myObjGrad(decVars, pack,dynamics, pathObj, bndObj, gradInfo)
%
% This function unpacks the decision variables, sends them to the
% user-defined objective functions, and then returns the final cost
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   bndObj = user-defined boundary objective function
%   gradInfo =
%
% OUTPUTS:
%   cost = scalar cost for this set of decision variables
%   dcost = gradient of cost
%     NOTE: gradients are only available for pathCost that depends only on
%     input parameters not states.
%
%

% All of the real work happens inside this function:
[t,x,~,~,pathCost,dxdalpha,dJdalpha] = simSysGrad(decVars, pack, dynamics, pathObj, gradInfo); %#ok<ASGLU>
% dxdalpha is included in outputs to make sure subsequent calls to
% simulateSystem without change a to decVars have access to the correct value
% of dxdalpha - see simulateSystem in which dxdalpha is not calculated unless
% nargout > 5

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

cost = pathCost + bndCost;

% calculate gradient of cost function
if nargout > 1
    
    nState = pack.nState;
    nControl = pack.nControl;
    nSegment = pack.nSegment;
    nSubStep = pack.nSubStep;
    nDecVar = 2+nState*(1+nSegment)+nControl*(1+nSegment*nSubStep*2);
    
    % allocate gradient of cost
    dcost_pth = zeros(nDecVar,1);
    dcost_bnd = zeros(nDecVar,1);
    
    % gradient assocated with bound objective
    if ~isempty(bndObj)
        
        % bound costs and gradients w.r.t. t0, x0, tF, xF
        [~, d_bnd] = bndObj(t0,x0,tF,xF);
        
        % gradients of t0, x0, tF, xF w.r.t. decision parameters (labeled alpha)
        dt0_dalpha = zeros(1,nDecVar);
        dt0_dalpha(1) = 1; % t0 is always the first decVar
        %
        dx0_dalpha = zeros(nState,nDecVar);
        dx0_dalpha(1:nState,gradInfo.xIdx(:,end)) = eye(nState);
        %
        dtF_dalpha = zeros(1,nDecVar);
        dtF_dalpha(2) = 1; % tF is always the second decVar
        %
        dxF_dalpha = zeros(nState,nDecVar);
        dxF_dalpha(1:nState,gradInfo.xIdx(:,end)) = eye(nState);
        
        % gradient of bound cost
        dcost_bnd(:) = [dt0_dalpha; dx0_dalpha; dtF_dalpha; dxF_dalpha]' * d_bnd';
    end
    
    % gradient assocated with path objective
    if ~isempty(pathObj)
        
        dcost_pth = dJdalpha';
        
    end
    
    dcost = dcost_pth + dcost_bnd;
    
end

end

function [c, ceq, dc, dceq] = myCstGrad(decVars, pack, dynamics, pathObj, pathCst, bndCst, gradInfo)
%
% This function unpacks the decision variables, computes the defects along
% the trajectory, and then evaluates the user-defined constraint functions.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   pathCst = user-defined path-constraint function
%   bndCst = user-defined boundary constraint function
%   gradInfo =
%
% OUTPUTS:
%   c = non-linear inequality constraint
%   ceq = non-linear equatlity cosntraint
%   dc = gradient of c w.r.t. decVars
%   dceq = gradient of ceq w.r.t. decVars
%
% NOTE:
%   - path constraints are  satisfied at the start and end of each sub-step
%


[t,x,u,defects,pathcost,dxdalpha] = simSysGrad(decVars, pack, dynamics, pathObj, gradInfo); %#ok<ASGLU>

%%%% Call user-defined constraints and pack up:
if nargout <= 2
    [c, ceq] = collectConstraints(t,x,u,...
        defects,...
        pathCst, bndCst);
    
else
    
    [c, ceq, dc, dceq] = collectConstraintsGrad(t,x,u,...
        defects,...
        pathCst, bndCst, pack, gradInfo, dxdalpha);
    
end

end


function [c, ceq, dc, dceq] = collectConstraintsGrad(t,x,u,defects, pathCst, bndCst, pack, gradInfo, dxdalpha)
% [c, ceq, dc, dceq] = collectConstraints(t,x,u,defects, pathCst, bndCst, pack, gradInfo, dxdalpha)
%
% TrajOpt utility function.
%
% Collects the defects, calls user-defined constraints, and then packs
% everything up into a form that is good for fmincon.
%
% INPUTS:
%   t = time vector (time at each substep) nTime = 1+nSegment*nSubStep
%   x = state matrix (states at each time in t)
%   u = control matrix (control at each time in t)
%   defects = defects matrix
%   pathCst = user-defined path constraint function
%   bndCst = user-defined boundary constraint function
%   pack =
%   gradInfo =
%   dxdalpha = partial derivative of state at each substep w.r.t. decVars
%
% OUTPUTS:
%   c = inequality constraint for fmincon
%   ceq = equality constraint for fmincon
%   dc = gradient of c w.r.t. decVars
%   dceq = gradient of ceq w.r.t. decVars
%

% problem dimensions
nState = pack.nState;
nControl = pack.nControl;
nSegment = pack.nSegment;
nSubStep = pack.nSubStep;
nDecVar = 2+nState*(1+nSegment)+nControl*(1+nSegment*nSubStep*2);

%%%% defect constraints
ceq_dyn = reshape(defects,numel(defects),1);

dceq_dyn = zeros(nDecVar,length(ceq_dyn));
Inx = eye(nState);
for j = 1:nSegment
    rows = (j-1)*nState+(1:nState);
    cols = rows;
    dceq_dyn(:,cols) = dxdalpha{j}(:,:,end)';  % gradient w.r.t. to x_i(+)
    dceq_dyn(2+nState+rows,cols) = -Inx; % gradient w.r.t. to x_i
end


%%%% Compute the user-defined constraints:

%%%% path constraints
if isempty(pathCst)
    c_path = [];
    ceq_path = [];
    dc_path = [];
    dceq_path = [];
else
    [c_pathRaw, ceq_pathRaw, c_pathGradRaw, ceq_pathGradRaw] = pathCst(t,x,u);
    c_path = reshape(c_pathRaw,numel(c_pathRaw),1);
    ceq_path = reshape(ceq_pathRaw,numel(ceq_pathRaw),1);
    
    dc_path = zeros(nDecVar,length(c_path));
    dceq_path = zeros(nDecVar,length(ceq_path));
    
    % dt/dalpha : gradient of time w.r.t. decVars
    dt_dalpha = zeros(1,nDecVar);
    nTime = 1+nSegment*nSubStep;
    n_time = 0:nTime-1;
    
    % gradients of path constraints
    nc = size(c_pathRaw,1); % number path constraints at each time
    nceq = size(ceq_pathRaw,1);
    for j = 1:(nSegment+1)
        for i = 1:nSubStep
            
            % d(t[n])/dalpha
            n_time0 = n_time((j-1)*nSubStep+i);
            dt_dalpha(1) = (1 - n_time0/(nTime-1));
            dt_dalpha(2) = (n_time0/(nTime-1));
            
            %
            if j < nSegment+1
                dxi_dalpha = dxdalpha{j}(:,:,i);
            else
                dxi_dalpha = zeros(nState,nDecVar);
                cols = gradInfo.xIdx(:,j);
                dxi_dalpha(:,cols) = eye(nState);
            end
            
            %
            dui_dalpha = zeros(nControl,nDecVar);
            cols = gradInfo.indu(:,(j-1)*nSubStep+i);
            dui_dalpha(:,cols) = eye(nControl);
            
            % inequality path constraints
            if nc > 0
                cols = (1:nc) + nc*((j-1)*nSubStep+i-1);
                dc_path(:,cols) = [dt_dalpha; dxi_dalpha; dui_dalpha]' * c_pathGradRaw(:,:,nSubStep*(j-1)+i)';
            end
            
            % equality path constraints
            if nceq > 0
                cols = (1:nceq) + nceq*((j-1)*nSubStep+i-1);
                dceq_path(:,cols) = [dt_dalpha; dxi_dalpha; dui_dalpha]' * ceq_pathGradRaw(:,:,nSubStep*(j-1)+i)';
            end
            
            % no need to continue with inner loop.
            if j == nSegment+1
                break;
            end
        end
    end
    
end

%%%% bound constraints
if isempty(bndCst)
    c_bnd = [];
    ceq_bnd = [];
    dc_bnd = [];
    dceq_bnd = [];
    
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    
    % bound constraints and gradients w.r.t. t0, x0, tF, xF
    [c_bnd, ceq_bnd, d_bnd, deq_bnd] = bndCst(t0,x0,tF,xF);
    
    % gradients of t0, x0, tF, xF w.r.t. decision parameters (labeled alpha)
    dt0_dalpha = zeros(1,nDecVar);
    dt0_dalpha(1) = 1; % t0 is always the first decVar
    %
    dx0_dalpha = zeros(nState,nDecVar);
    cols = 2+(1:nState);
    dx0_dalpha(1:nState,cols) = eye(nState);
    %
    dtF_dalpha = zeros(1,nDecVar);
    dtF_dalpha(2) = 1; % tF is always the second decVar
    %
    dxF_dalpha = zeros(nState,nDecVar);
    cols = (1:nState) + 2 + nSegment*nState;
    dxF_dalpha(1:nState,cols) = eye(nState);
    
    
    % inequality bound constraints
    dc_bnd = [dt0_dalpha; dx0_dalpha; dtF_dalpha; dxF_dalpha]' * d_bnd';
    
    % equality bound constraints
    dceq_bnd = [dt0_dalpha; dx0_dalpha; dtF_dalpha; dxF_dalpha]' * deq_bnd';
    
end

%%%% Pack everything up:
c = [c_path;c_bnd];
ceq = [ceq_dyn; ceq_path; ceq_bnd];

dc = [dc_path, dc_bnd];
dceq = [dceq_dyn, dceq_path, dceq_bnd];


end



%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [t,x,u,defects,pathCost,dxdalpha,dJdalpha] = simSysGrad(decVars, pack, dynFun, pathObj, gradInfo)
%
% This function does the real work of the transcription method. It
% simulates the system forward in time across each segment of the
% trajectory, computes the integral of the cost function, and then matches
% up the defects between the end of each segment and the start of the next.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%
% OUTPUTS:
%   t = [1 x nGrid] = time vector for the edges of the sub-step grid
%   x = [nState x nGrid] = state vector
%   u = [nControl x nGrid] = control vector
%   defects = [nState x nSegment] = defect matrix
%   pathCost = scalar cost for the path integral
%
% NOTES:
%   - nGrid = nSegment*nSubStep+1
%   - This function is usually called twice for each combination of
%   decision variables: once by the objective function and once by the
%   constraint function. To keep the code fast I cache the old values and
%   only recompute when the inputs change.
%


%%%% CODE OPTIMIZATION %%%%
%
% Prevents the same exact code from being called twice by caching the
% solution and reusing it when appropriate.
%
global RUNGE_KUTTA_t RUNGE_KUTTA_x RUNGE_KUTTA_u
global RUNGE_KUTTA_defects RUNGE_KUTTA_pathCost
global RUNGE_KUTTA_decVars RUNGE_KUTTA_dzdalpha
%
usePreviousValues = false;
if ~isempty(RUNGE_KUTTA_decVars)
    if length(RUNGE_KUTTA_decVars) == length(decVars)
        if ~any(RUNGE_KUTTA_decVars ~= decVars)
            usePreviousValues = true;
        end
    end
end
%
if usePreviousValues
    t = RUNGE_KUTTA_t;
    x = RUNGE_KUTTA_x;
    u = RUNGE_KUTTA_u;
    defects = RUNGE_KUTTA_defects;
    pathCost = RUNGE_KUTTA_pathCost;
    dxdalpha = RUNGE_KUTTA_dzdalpha;
else
    %
    %
    %%%% END CODE OPTIMIZATION %%%%
    
    
    [tSpan, state, control] = unPackDecVar(decVars,pack);
    
    nState = pack.nState;
    nControl = pack.nControl;
    nSegment = pack.nSegment;
    nSubStep = pack.nSubStep;
    
    % NOTES:
    %   The following bit of code is a bit confusing, mostly due to the
    %   need for vectorization to make things run at a reasonable speed in
    %   Matlab. Part of the confusion comes because the decision variables
    %   include the state at the beginning of each segment, but the control
    %   at the beginning and middle of each substep - thus there are more
    %   control grid-points than state grid points. The calculations are
    %   vectorized over segments, but not sub-steps, since the result of
    %   one sub-step is required for the next.
    
    % time, state, and control at the ends of each substep
    nTime = 1+nSegment*nSubStep;
    t = linspace(tSpan(1), tSpan(2), nTime);
    x = zeros(nState, nTime);
    u = control(:,1:2:end); % Control a the endpoints of each segment
    uMid = control(:,2:2:end);  %Control at the mid-points of each segment
    c = zeros(1, nTime-1);  %Integral cost for each segment
    dt = (t(end)-t(1))/(nTime-1);
    
    idx = 1:nSubStep:(nTime-1);   %Indicies for the start of each segment
    x(:,[idx,end]) = state;   %Fill in the states that we already know
    
    % VARIABLES for analytic gradient evaluations.
    % size of decicion parameters (2 for time), nstate*(nSegment+1), ...
    % dxdalpha = partial derivative of state w.r.t. decVars (alpha)
    nalpha = 2 + nState*(1+nSegment) + nControl*(1+2*nSubStep*nSegment);
    dxdalpha = cell(1,nSegment);
    for i = 1:nSegment
        dxdalpha{i} = zeros(nState,nalpha,nSubStep+1);
        cols = 2+(i-1)*nState+(1:nState);
        dxdalpha{i}(:,cols,1) = eye(nState);
    end
    dTdalpha = zeros(1,nalpha); dTdalpha(1:2) = [-1,1];
    dt_dalpha = zeros(1,nalpha);
    n_time = 0:nTime-1;
    
    % gradient of path cost
    dJdalpha = zeros(1,nalpha);
    
    for iSubStep = 1:nSubStep
        % March forward Runge-Kutta step
        
        t0 = t(idx);
        x0 = x(:,idx);
        
        
        
        %------------------------------------------
        % Code for calculating dxdalpha (partial derivative of state w.r.t.
        % the descision parameters): dxdalpha = nstate x nalpha
        % assume nargout <=5 when using finite difference calculation for
        % gradients in which case dxdalpha is unnecessary.
        
        % Gradient of time w.r.t. decVars
        % ------------------------------------------------------------
        % dt = (tF-t0)/(nTime-1)
        % t = t0 + n*dt
        % t = t0 + n*(tF-t0)/(nTime-1)
        % t = t0*(1-n/(nTime-1)) + tF*(n/(nTime-1))
        %
        % alpha = [t0, tF, x0, x1, ..., xN, u0, uM0, u1, ..., uN]
        % dt/dalpha = [1 - n/(nTime-1), n/(nTime-1), 0, 0, ... 0]
        % ------------------------------------------------------------
        
        n_time0 = n_time(idx);
        
        [k0, dk0] = combinedDynGrad(t0,        x0,                         u(:,idx), dynFun,pathObj);
        [k1, dk1] = combinedDynGrad(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid(:,idx), dynFun,pathObj);
        [k2, dk2] = combinedDynGrad(t0+0.5*dt, x0 + 0.5*dt*k1(1:nState,:), uMid(:,idx), dynFun,pathObj);
        [k3, dk3] = combinedDynGrad(t0+dt,     x0 +     dt*k2(1:nState,:), u(:,idx+1), dynFun,pathObj);
        z = (dt/6)*(k0 + 2*k1 + 2*k2 + k3);  %Change over the sub-step
        
        for j = 1:nSegment
            
            % d(t[n])/dalpha
            dt_dalpha(1) = (1 - n_time0(j)/(nTime-1));
            dt_dalpha(2) = (n_time0(j)/(nTime-1));
            
            % du[n]/dalpha
            du_dalpha = zeros(nControl,nalpha);
            du_dalpha(:,gradInfo.indu(:,idx(j))) = eye(nControl);
            
            % duMid[n]/dalpha
            duMid_dalpha = zeros(nControl,nalpha);
            duMid_dalpha(:,gradInfo.indumid(:,idx(j))) = eye(nControl);
            
            % du[n+1]/dalpha
            du1_dalpha = zeros(nControl,nalpha);
            du1_dalpha(:,gradInfo.indu(:,idx(j)+1)) = eye(nControl);
            
            % dk0/dalpha
            dk0da = dk0(:,:,j) * [dt_dalpha; dxdalpha{j}(:,:,iSubStep); du_dalpha];
            
            % dk1/dalpha
            dk1da = dk1(:,:,j) * [dt_dalpha + 0.5/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + 0.5*dt*dk0da(1:nState,:) + 0.5/(nTime-1)*k0(1:nState,j)*dTdalpha; duMid_dalpha];
            
            % dk2/dalpha
            dk2da = dk2(:,:,j) * [dt_dalpha + 0.5/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + 0.5*dt*dk1da(1:nState,:) + 0.5/(nTime-1)*k1(1:nState,j)*dTdalpha; duMid_dalpha];
            
            % dk3/dalpha
            dk3da = dk3(:,:,j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + dt*dk2da(1:nState,:) + 1/(nTime-1)*k2(1:nState,j)*dTdalpha; du1_dalpha];
            
            dz = (dt/6)*(dk0da + 2*dk1da + 2*dk2da + dk3da)...
                + 1/(6*(nTime-1))*(k0(:,j)+2*k1(:,j)+2*k2(:,j)+k3(:,j))*dTdalpha;
            
            % update dxdalpha
            dxdalpha{j}(:,:,iSubStep+1) = dxdalpha{j}(:,:,iSubStep) + dz(1:nState,:);
            
            % update dJdalpha
            dJdalpha  = dJdalpha + dz(nState+1,:);
        end
        
        
        xNext = x0 + z(1:nState,:);  %Next state
        c(idx) = z(end,:);  %Integral of the cost function over this step
        
        if iSubStep == nSubStep %We've reached the end of the interval
            % Compute the defect vector:
            defects = xNext - x(:,idx+1);
        else
            % Store the state for next step in time
            idx = idx+1;   %  <-- This is important!!
            x(:,idx) = xNext;
        end
        
    end
    
    pathCost = sum(c);  %Sum up the integral cost over each segment
    
    %%%% Cache results to use on the next call to this function.
    RUNGE_KUTTA_t = t;
    RUNGE_KUTTA_x = x;
    RUNGE_KUTTA_u = u;
    RUNGE_KUTTA_defects = defects;
    RUNGE_KUTTA_pathCost = pathCost;
    
end

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [dz, J] = combinedDynGrad(t,x,u,dynFun,pathObj)
% [dz, dJ] = combinedDynGrad(t,x,u,dynFun,pathObj)
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
%   dJ = [JAC(dynamics), JAC(objective)] = combined jacobian of dynamics
%   and objective w.r.t. (t,x,u)



nState = size(x,1);
nControl = size(u,1);

[dx,Jx] = dynFun(t,x,u);
if isempty(pathObj)
    dc = zeros(size(t));
    Jc = zeros(1,1+nState+nControl,length(t));
else
    [dc,Jc] = pathObj(t,x,u);
    Jc = reshape(Jc,1,1+nState+nControl,length(t));
end

dz = [dx;dc];

J = cat(1,Jx,Jc);


end
