function soln = DEV_dirColShooting(problem)
% soln = directCollocation(problem)
%
% TrajOpt utility function
%
% This function is designed to be called by either "trapezoid" or
% "hermiteSimpson". It actually calls FMINCON to solve the trajectory
% optimization problem. 
%
%

%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

nGrid = length(F.weights);

flagGradObj = strcmp(Opt.nlpOpt.GradObj,'on');
flagGradCst = strcmp(Opt.nlpOpt.GradConstr,'on');

% Print out notice about analytic gradients
if Opt.verbose > 0
    if flagGradObj
        fprintf('      - using analytic gradients of objective function\n');
    end
    if flagGradCst
        fprintf('      - using analytic gradients of constraint function\n');
    end
    fprintf('\n');
end


switch Opt.method
    case 'trapezoid'
        % number of shooting segments
        nShootSegment = Opt.trapezoid.nShootSegment;
        % number of splines per shooting Segment
        nSplineSegment = floor((nGrid-1)/nShootSegment);
        % number of grid points per shooting segment
        nGridSegment = 1+nSplineSegment; 
        % State indices corresponding to the end of a shooting segment
        idx_ShootEnd = nGridSegment : nGridSegment : nShootSegment*nGridSegment;
        %idx_ShootEnd = nGridSegment : nGridSegment : nGrid+nShootSegment;
        % Main Trajectory indices
        idx_Traj = 1:nGrid+nShootSegment;
        idx_Traj(idx_ShootEnd) = [];
        
        % the last segment has more points than other segments
        if (idx_Traj(end) - (idx_ShootEnd(end)+1) > nSplineSegment)
            warning(['The number of splines in the last shooting ',...
                'segment is greater than the other segments. May be good to reduce the number of shooting segments.'])
            fprintf('\n')
        end
        
    case 'hermiteSimpson'
        % number of shooting segments
        nShootSegment = Opt.hermiteSimpson.nShootSegment;
        % number of splines per shooting Segment
        nSplineSegment = floor((Opt.hermiteSimpson.nSegment)/nShootSegment);
        % number of grid points per shooting segment
        nGridSegment = 1+2*nSplineSegment; 
        % State indices corresponding to the end of a shooting segment
        idx_ShootEnd = nGridSegment : nGridSegment : nShootSegment*nGridSegment;
        %idx_ShootEnd = nGridSegment : nGridSegment : nGrid+nShootSegment;
        % Main Trajectory indices
        idx_Traj = 1:nGrid+nShootSegment;
        idx_Traj(idx_ShootEnd) = [];
        
        % the last segment has more points than other segments
        if (idx_Traj(end) - (idx_ShootEnd(end)+1) > nSplineSegment)
            warning(['The number of splines in the last shooting ',...
                'segment is greater than the other segments. May be good to reduce the number of shooting segments.'])
            fprintf('\n')
        end
    
    otherwise
        error('undefined dirCol Shooting Method')
end


if Opt.verbose > 0
    fprintf('  -> Number of shooting segments = %d\n',nShootSegment);
    fprintf('  -> Number of splines per segment = %d\n',nSplineSegment);
    fprintf('\n');
end

% Interpolate the guess at the grid-points for transcription:
guess.tSpan = G.time([1,end]);

guess.time = linspace(guess.tSpan(1), guess.tSpan(2), nGrid);
% guess of state: Set x_minus = x_plus for each shooting segment
guess.state = zeros(size(G.state,1),nGrid+nShootSegment);
guess.state(:,idx_Traj) = interp1(G.time', G.state', guess.time')';
guess.state(:,idx_ShootEnd) = guess.state(:,idx_ShootEnd+1);
% guess of control: u_minus = u_plus
guess.control = zeros(size(G.control,1),nGrid+nShootSegment);
guess.control(:,idx_Traj) = interp1(G.time', G.control', guess.time')';
guess.control(:,idx_ShootEnd) = guess.control(:,idx_ShootEnd+1);

[zGuess, pack] = packDecVar(guess.time, guess.state, guess.control);

%%%% DirCol Shooting Info
% grid points associated with the main trajectory
pack.idx_Traj = idx_Traj;
% grid points associated with the shooting segment endpoints
pack.idx_ShootEnd = idx_ShootEnd;
% number of shooting segments
pack.nShootSegment = nShootSegment;


if flagGradCst || flagGradObj
    gradInfo = grad_computeInfo(pack);
end


% Unpack all bounds:
tLow = linspace(B.initialTime.low, B.finalTime.low, nGrid);
xLow = [B.initialState.low, B.state.low*ones(1,nGrid+nShootSegment-2), B.finalState.low];
uLow = B.control.low*ones(1,nGrid+nShootSegment);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = linspace(B.initialTime.upp, B.finalTime.upp, nGrid);
xUpp = [B.initialState.upp, B.state.upp*ones(1,nGrid+nShootSegment-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nGrid+nShootSegment);
zUpp = packDecVar(tUpp,xUpp,uUpp);

%%%% Set up problem for fmincon:
if flagGradObj
    P.objective = @(z)( ...
        myObjGrad(z, pack, F.pathObj, F.bndObj, F.weights, gradInfo) );   %Analytic gradients
else
    P.objective = @(z)( ...
        myObjective(z, pack, F.pathObj, F.bndObj, F.weights) );   %Numerical gradients
end
if flagGradCst
    P.nonlcon = @(z)( ...
        myCstGrad(z, pack, F.dynamics, F.pathCst, F.bndCst, F.defectCst, gradInfo) ); %Analytic gradients
else
    P.nonlcon = @(z)( ...
        myConstraint(z, pack, F.dynamics, F.pathCst, F.bndCst, F.defectCst) ); %Numerical gradients
end


% Plot Defect Constraint sparsity
if strcmp(Opt.(Opt.method).PlotDefectGrad,'on')
    [~,~,~,dceq] = myCstGrad(zGuess, pack, F.dynamics, [], [], F.defectCst, gradInfo);
    figure(100),clf
    spy(dceq')
end


% Check analytic gradients with DERIVEST package
if strcmp(Opt.(Opt.method).AdaptiveDerivativeCheck,'on')
    if exist('jacobianest','file')
        runGradientCheck(zGuess, pack,F.dynamics, F.pathObj, F.bndObj, F.weights, F.pathCst, F.bndCst, F.defectCst, gradInfo);
        Opt.nlpOpt.DerivativeCheck = [];  %Disable built-in derivative check
    else
        Opt.(Opt.method).AdaptiveDerivativeCheck = 'cannot find jacobianest.m';
        disp('Warning: the derivest package is not on search path.');
        disp(' --> Using fmincon''s built-in derivative checks.');
    end
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

% number of shoooting segments
nShootSegment = size(x,2) - nTime;

tSpan = [t(1); t(end)];
xCol = reshape(x, nState*(nTime+nShootSegment), 1);
uCol = reshape(u, nControl*(nTime+nShootSegment), 1);

% control and state indices in the decision variable vector
indz = reshape(2+(1:numel(u)+numel(x)),nState+nControl,nTime+nShootSegment);

% index of time, state, control variables in the decVar vector
tIdx = 1:2;
xIdx = indz(1:nState,:);
uIdx = indz(nState+(1:nControl),:);

% decision variables
% variables are indexed so that the defects gradients appear as a banded
% matrix
z = zeros(2+numel(indz),1);
z(tIdx(:),1) = tSpan;
z(xIdx(:),1) = xCol;
z(uIdx(:),1) = uCol;

pack.nTime = nTime;
pack.nState = nState;
pack.nControl = nControl;
pack.nDecVar = 2 + nState*(nTime+nShootSegment) + nControl*(nTime+nShootSegment);
pack.tIdx = tIdx;
pack.xIdx = xIdx;
pack.uIdx = uIdx;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [t,x,u,t_all,x_all,u_all] = unPackDecVar(z,pack)
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

idx_Traj = pack.idx_Traj;
idx_ShootEnd = pack.idx_ShootEnd;
nShootSegment = pack.nShootSegment;

nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;

t = linspace(z(1),z(2),nTime);

% tmp_u = control at trajectory points and shooting segment end points
tmp_u = reshape(z(pack.uIdx),nControl,nTime+nShootSegment);
% Control values at trajectory points, not shooting segment end points
u = tmp_u(:,idx_Traj);

% tmp_x = state at trajectory points and shooting segment end points
tmp_x = reshape(z(pack.xIdx),nState,nTime+nShootSegment);
% State values at trajectory points, not shooting segment end points
x = tmp_x(:,idx_Traj);

% Values of t,x,u at the end of shooting segments.
if nargout > 3
  
    % state and control at all points (trajectory + shooting segment end
    % points)
    u_all = tmp_u;
    x_all = tmp_x;

    % the time indices at the end of shooting segments
    idt_SegmentEnd = idx_ShootEnd - (0:nShootSegment-1);
    
    % create a time vector that includes times at the Shooting Endpoints
    t_all = zeros(1,nTime+nShootSegment);
    t_all(1,idx_Traj) = t;
    t_all(1,idx_ShootEnd) = t(1,idt_SegmentEnd);
    
end

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function cost = myObjective(z,pack,pathObj,bndObj,weights)
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

function [c, ceq] = myConstraint(z,pack,dynFun, pathCst, bndCst, defectCst)
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

[t,x,u,t_all,x_all,u_all] = unPackDecVar(z,pack);

%%%% Compute defects along the trajectory:
dt = (t(end)-t(1))/(length(t)-1);
f_all = dynFun(t_all,x_all,u_all);
defects = defectCst(pack,dt,x_all,u_all,f_all);

%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,defects, pathCst, bndCst);

end




%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%            Additional Sub-Functions for Gradients                 %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%

function [fail] = runGradientCheck(z_test, pack,dynamics, pathObj, bndObj, weights, pathCst, bndCst, defectCst, gradInfo)
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
[~, dcost] = myObjGrad(z_test, pack, pathObj, bndObj, weights, gradInfo);

% check gradient with derivest package
deriv = gradest(@(z) myObjGrad(z, pack, pathObj, bndObj, weights, gradInfo),z_test);

% print largest difference in numerical and analytic gradients
fprintf('\n%s\n','Objective function derivatives:')
fprintf('%s\n','Maximum relative difference between user-supplied')
fprintf('%s %1.5e \n','and finite-difference derivatives = ',max(abs(dcost-deriv')))
if any(abs(dcost-deriv') > GradientCheckTol)
    error('Objective gradient did not pass')
end

% analytic nonlinear constraints
[c, ceq,dc, dceq] = myCstGrad(z_test, pack, dynamics, pathCst, bndCst, defectCst, gradInfo);

% check nonlinear inequality constraints with 'jacobianest'
if ~isempty(c)
    jac = jacobianest(@(z) myCstGrad(z, pack, dynamics, pathCst, bndCst, defectCst, gradInfo),z_test);
    
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
    jac = jacobianest(@(z) myCstGradCheckEq(z, pack, dynamics, pathCst, bndCst, defectCst, gradInfo),z_test);
    
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

%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%

function ceq = myCstGradCheckEq(decVars, pack, dynamics, pathCst, bndCst, defectCst, gradInfo)
% This function is necessary for runGradientCheck function
% return only equality constraint (ceq) for use with jacobest.m

[~, ceq] = myCstGrad(decVars, pack, dynamics, pathCst, bndCst, defectCst, gradInfo);

end


%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%


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


idx_Traj = pack.idx_Traj;
idx_ShootEnd = pack.idx_ShootEnd;
nShootSegment = pack.nShootSegment;

nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;
nDecVar = pack.nDecVar;

gradInfo.nDecVar = nDecVar;
%zIdx = 1:nDecVar;
%[tIdx, xIdx, uIdx, tIdxAll, xIdxAll, uIdxAll] = unPackDecVar(zIdx,pack);
tIdx = pack.tIdx;
xIdx = pack.xIdx(:,pack.idx_Traj);
uIdx = pack.uIdx(:,pack.idx_Traj);
xIdxAll = pack.xIdx;
uIdxAll = pack.uIdx;

% indices at the trajectory points
gradInfo.tIdx = tIdx([1,end]);
gradInfo.xuIdx = [xIdx;uIdx];

% indices at trajectory points and Shooting Segment Endpoints
gradInfo.xuIdxAll = [xIdxAll;uIdxAll];

%%%% Compute gradients of time:
% alpha = (0..N-1)/(N-1)
% t = alpha*tUpp + (1-alpha)*tLow
alpha = (0:(nTime-1))/(nTime-1);
gradInfo.alpha = [1-alpha; alpha];

% create an alpha vector that includes alpha at the Shooting Endpoints
alphaAll = zeros(1,nTime+nShootSegment);
alphaAll(1,idx_Traj) = alpha;
idt_SegmentEnd = idx_ShootEnd - (0:nShootSegment-1);
alphaAll(1,idx_ShootEnd) = alpha(1,idt_SegmentEnd);
gradInfo.alphaAll = [1-alphaAll;alphaAll];

if (gradInfo.tIdx(1)~=1 || gradInfo.tIdx(end)~=2)
    error('The first two decision variables must be the initial and final time')
end
gradInfo.dtGrad = [-1; 1]/(nTime-1);

%%%% Compute gradients of state
gradInfo.xGrad = zeros(nState,nTime,nDecVar);
for iTime=1:nTime
    for iState=1:nState
        gradInfo.xGrad(iState,iTime,xIdx(iState,iTime)) = 1;
    end
end

%%%% Compute gradients of state including the decVars of states at shooting
%%%% Segment Endpoints 
gradInfo.xGradAll = zeros(nState,nTime,nDecVar);
for iTime=1:nTime+nShootSegment
    for iState=1:nState
        gradInfo.xGradAll(iState,iTime,xIdxAll(iState,iTime)) = 1;
    end
end

gradInfo.uGradAll = zeros(nControl,nTime,nDecVar);
for iTime=1:nTime+nShootSegment
    for iControl=1:nControl
        gradInfo.uGradAll(iControl,iTime,uIdxAll(iControl,iTime)) = 1;
    end
end

%%%% For unpacking the boundary constraints and objective:
gradInfo.bndIdxMap = [tIdx(1); xIdx(:,1); tIdx(end); xIdx(:,end)];


end


%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%


function [c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,defects, defectsGrad, pathCst, bndCst, gradInfo)
% [c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,defects, defectsGrad, pathCst, bndCst, gradInfo)
%
% TrajOpt utility function.
%
% Collects the defects, calls user-defined constraints, and then packs
% everything up into a form that is good for fmincon. Additionally, it
% reshapes and packs up the gradients of these constraints.
%
% INPUTS:
%   t = time vector
%   x = state matrix
%   u = control matrix
%   defects = defects matrix
%   pathCst = user-defined path constraint function
%   bndCst = user-defined boundary constraint function
%
% OUTPUTS:
%   c = inequality constraint for fmincon
%   ceq = equality constraint for fmincon
%

ceq_dyn = reshape(defects,numel(defects),1);
if size(defectsGrad,3) > 1
    ceq_dynGrad = grad_flattenPathCst(defectsGrad);
else
    ceq_dynGrad = defectsGrad;
end

%%%% Compute the user-defined constraints:
if isempty(pathCst)
    c_path = [];
    ceq_path = [];
    c_pathGrad = [];
    ceq_pathGrad = [];
else
    [c_pathRaw, ceq_pathRaw, c_pathGradRaw, ceq_pathGradRaw] = pathCst(t,x,u);
    c_path = reshape(c_pathRaw,numel(c_pathRaw),1);
    ceq_path = reshape(ceq_pathRaw,numel(ceq_pathRaw),1);
    c_pathGrad = grad_flattenPathCst(grad_reshapeContinuous(c_pathGradRaw,gradInfo));
    ceq_pathGrad = grad_flattenPathCst(grad_reshapeContinuous(ceq_pathGradRaw,gradInfo));
end
if isempty(bndCst)
    c_bnd = [];
    ceq_bnd = [];
    c_bndGrad = [];
    ceq_bndGrad = [];
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    [c_bnd, ceq_bnd, c_bndGradRaw, ceq_bndGradRaw] = bndCst(t0,x0,tF,xF);
    c_bndGrad = grad_reshapeBoundary(c_bndGradRaw,gradInfo);
    ceq_bndGrad = grad_reshapeBoundary(ceq_bndGradRaw,gradInfo);
end

%%%% Pack everything up:
c = [c_path;c_bnd];
ceq = [ceq_dyn; ceq_path; ceq_bnd];

cGrad = [c_pathGrad;c_bndGrad]';
ceqGrad = [ceq_dynGrad; ceq_pathGrad; ceq_bndGrad]';


end



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



function C = grad_flattenPathCst(CC)
%
% This function takes a path constraint and reshapes the first two
% dimensions so that it can be passed to fmincon
%
if isempty(CC)
    C = [];
else
    [n1,n2,n3] = size(CC);
    C = reshape(CC,n1*n2,n3);
end

end





%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



function CC = grad_reshapeBoundary(C,gradInfo)
%
% This function takes a boundary constraint or objective from the user
% and expands it to match the full set of decision variables
%

CC = zeros(size(C,1),gradInfo.nDecVar);
CC(:,gradInfo.bndIdxMap) = C;

end



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



function grad = grad_reshapeContinuous(gradRaw,gradInfo)
% grad = grad_reshapeContinuous(gradRaw,gradInfo)
%
% TrajOpt utility function.
%
% This function converts the raw gradients from the user function into
% gradients with respect to the decision variables.
%
% INPUTS:
%   stateRaw = [nOutput,nInput,nTime]
%
% OUTPUTS:
%   grad = [nOutput,nTime,nDecVar]
%

if isempty(gradRaw)
    grad = [];
else
    [nOutput, ~, nTime] = size(gradRaw);
    
    grad = zeros(nOutput,nTime,gradInfo.nDecVar);
    
    % First, loop through and deal with time.
    timeGrad = gradRaw(:,1,:); timeGrad = permute(timeGrad,[1,3,2]);
    for iOutput=1:nOutput
        A = ([1;1]*timeGrad(iOutput,:)).*gradInfo.alpha;
        grad(iOutput,:,gradInfo.tIdx) = permute(A,[3,2,1]);
    end
    
    % Now deal with state and control:
    for iOutput=1:nOutput
        for iTime=1:nTime
            B = gradRaw(iOutput,2:end,iTime);
            grad(iOutput,iTime,gradInfo.xuIdx(:,iTime)) = permute(B,[3,1,2]);
        end
    end
end

end


function grad = grad_reshapeContinuousAll(gradRaw,gradInfo)
% grad = grad_reshapeContinuous(gradRaw,gradInfo)
%
% TrajOpt utility function.
%
% This function converts the raw gradients from the user function into
% gradients with respect to the decision variables.
%
% INPUTS:
%   stateRaw = [nOutput,nInput,nTime]
%
% OUTPUTS:
%   grad = [nOutput,nTime,nDecVar]
%

if isempty(gradRaw)
    grad = [];
else
    [nOutput, ~, nTime] = size(gradRaw);
    
    grad = zeros(nOutput,nTime,gradInfo.nDecVar);
    
    % First, loop through and deal with time.
    timeGrad = gradRaw(:,1,:); timeGrad = permute(timeGrad,[1,3,2]);
    for iOutput=1:nOutput
        A = ([1;1]*timeGrad(iOutput,:)).*gradInfo.alphaAll;
        grad(iOutput,:,gradInfo.tIdx) = permute(A,[3,2,1]);
    end
    
    % Now deal with state and control:
    for iOutput=1:nOutput
        for iTime=1:nTime
            B = gradRaw(iOutput,2:end,iTime);
            grad(iOutput,iTime,gradInfo.xuIdxAll(:,iTime)) = permute(B,[3,1,2]);
        end
    end
end

end



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



function [cost, costGrad] = myObjGrad(z,pack,pathObj,bndObj,weights,gradInfo)
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

%Unpack the decision variables:
[t,x,u] = unPackDecVar(z,pack);

% Time step for integration:
dt = (t(end)-t(1))/(length(t)-1);
dtGrad = gradInfo.dtGrad;
nTime = length(t);
nDecVar = length(z);

% Compute the cost integral along the trajectory
if isempty(pathObj)
    integralCost = 0;
    integralCostGrad = zeros(1,nDecVar);
else
        
    % Objective function integrand and gradients:
    [obj, objGradRaw] = pathObj(t,x,u);
    nInput = size(objGradRaw,1);
    objGradRaw = reshape(objGradRaw,1,nInput,nTime); 
    objGrad = grad_reshapeContinuous(objGradRaw,gradInfo);
    
    % integral objective function
    unScaledIntegral = obj*weights;
    integralCost = dt*unScaledIntegral;
    
    % Gradient of integral objective function
    dtGradTerm = zeros(1,nDecVar);
    dtGradTerm(1) = dtGrad(1)*unScaledIntegral;
    dtGradTerm(2) = dtGrad(2)*unScaledIntegral;
    objGrad = reshape(objGrad,nTime,nDecVar);
    integralCostGrad = ...
        dtGradTerm + ...
        dt*sum(objGrad.*(weights*ones(1,nDecVar)),1);
end

% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
    bndCostGrad = zeros(1,nDecVar);
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    [bndCost, bndCostGradRaw] = bndObj(t0,x0,tF,xF);
    bndCostGrad = grad_reshapeBoundary(bndCostGradRaw,gradInfo);    
end

% Cost function
cost = bndCost + integralCost;

% Gradients
costGrad = bndCostGrad + integralCostGrad;

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq, cGrad, ceqGrad] = myCstGrad(z,pack,dynFun, pathCst, bndCst, defectCst, gradInfo)
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


%Unpack the decision variables:
[t,x,u,t_all,x_all,u_all] = unPackDecVar(z,pack);

% Time step for integration:
dt = (t(end)-t(1))/(length(t)-1);
dtGrad = gradInfo.dtGrad;

% Gradient of the state with respect to decision variables
% xGrad = gradInfo.xGrad;
xGradAll = gradInfo.xGradAll;
uGradAll = gradInfo.uGradAll;

%%%% Compute defects along the trajectory:
% [f, fGradRaw] = dynFun(t,x,u);
% fGrad = grad_reshapeContinuous(fGradRaw,gradInfo);

[f_all, f_allGradRaw] = dynFun(t_all,x_all,u_all);
fGradAll = grad_reshapeContinuousAll(f_allGradRaw,gradInfo);

%computeDefects(dt,x,f_all,idx_ShootEnd,dtGrad,xGrad,fGrad)
[defects, defectsGrad] = defectCst(pack,dt,x_all,u_all,f_all,...
    dtGrad, xGradAll, uGradAll, fGradAll);

% Compute gradients of the user-defined constraints and then pack up:
[c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,...
    defects, defectsGrad, pathCst, bndCst, gradInfo);

end

