function soln = dirColGrad(P, problem)
% soln = dirColGrad(P, problem)
%
% TrajOpt utility function - Direct Collocation with Gradients
%
% This function is core function that is called to run the transcription
% for both the "trapezoid" and the "hermiteSimpson" methods when they are
% running analytic gradients.
%
% 

F = problem.func;

if isempty(P.objective)
    P.objective = @(z)( ...
        grad_objective(z, pack, F.pathObj, F.bndObj, gradInfo, weights) );   %Analytic gradients
end

if isempty(P.constraint)
     P.nonlcon = @(z)( ...
        myCstGrad(z, pack, F.dynamics, F.pathCst, F.bndCst, F.defectCst, gradInfo) ); %Analytic gradients
end


%%%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal,exitFlag,output] = fmincon(P);
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack);
nlpTime = toc;

%%%% Store the results:

soln.grid.time = tSoln;
soln.grid.state = xSoln;
soln.grid.control = uSoln;

soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct


end



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%


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



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



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


nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;
nDecVar = 2 + nState*nTime + nControl*nTime;

zIdx = 1:nDecVar;
gradInfo.nDecVar = nDecVar;
[tIdx, xIdx, uIdx] = unPackDecVar(zIdx,pack);
gradInfo.tIdx = tIdx([1,end]);
gradInfo.xuIdx = [xIdx;uIdx];

%%%% Compute gradients of time:
% alpha = (0..N-1)/(N-1)
% t = alpha*tUpp + (1-alpha)*tLow
alpha = (0:(nTime-1))/(nTime-1);
gradInfo.alpha = [1-alpha; alpha];

%%%% Compute gradients of state
gradInfo.xGrad = zeros(nState,nTime,nDecVar);
for iTime=1:nTime
    for iState=1:nState
        gradInfo.xGrad(iState,iTime,xIdx(iState,iTime)) = 1;
    end
end

%%%% For unpacking the boundary constraints and objective:
gradInfo.bndIdxMap = [tIdx(1); xIdx(:,1); tIdx(end); xIdx(:,end)];


end



%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%




function [dt, dtGrad] = grad_timeStep(t,gradInfo)
%
% TrajOpt utility function
%
% Computes the time step and its gradient
%
% dt = [1,1]
% dtGrad = [1,nz]
%

nTime = length(t);

dt = (t(end)-t(1))/(nTime-1);
dtGrad = zeros(1,gradInfo.nDecVar);

dtGrad(1,gradInfo.tIdx(1)) = -1/(nTime-1);
dtGrad(1,gradInfo.tIdx(2)) = 1/(nTime-1);

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
ceq_dynGrad = grad_flattenPathCst(defectsGrad);

%%%% Compute the user-defined constraints:
if isempty(pathCst)
    c_path = [];
    ceq_path = [];
    c_pathGrad = [];
    ceq_pathGrad = [];
else
    [c_path, ceq_path, c_pathGradRaw, ceq_pathGradRaw] = pathCst(t,x,u);
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






%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%


function [cost, costGrad] = grad_objective(z,pack,pathObj,bndObj,gradInfo,weights)
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
[dt, dtGrad] = grad_timeStep(t, gradInfo);
nTime = length(t);
nState = size(x,1);
nControl = size(u,1);
nDecVar = length(z);

% Compute the cost integral along the trajectory
if isempty(pathObj)
    integralCost = 0;
    integralCostGrad = zeros(nState+nControl,1);
else
    
    % Objective function integrand and gradients:
    [obj, objGradRaw] = pathObj(t,x,u);
    nInput = size(objGradRaw,1);
    objGradRaw = reshape(objGradRaw,1,nInput,nTime);
    objGrad = grad_reshapeContinuous(objGradRaw,gradInfo);
    
    % Integration by quadrature:
    integralCost = dt*obj*weights;  % Integration
    
    % Gradient of integral objective function
    objGrad = reshape(objGrad,nTime,nDecVar);
    integralCostGrad = ...
        dtGrad*(obj*weights) + ...
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







%%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%%



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
[t,x,u] = unPackDecVar(z,pack);

% Time step for integration:
[dt, dtGrad] = grad_timeStep(t,gradInfo);

%%%% Compute defects along the trajectory:
[f, fGradRaw] = dynFun(t,x,u);
fGrad = grad_reshapeContinuous(fGradRaw,gradInfo);

[defects, defectsGrad] = defectCst(dt,x,f,dtGrad,xGrad,fGrad,gradInfo);

% Compute gradients of the user-defined constraints and then pack up:
[c, ceq, cGrad, ceqGrad] = grad_collectConstraints(t,x,u,...
    defects, defectsGrad, pathCst, bndCst, gradInfo);

end


