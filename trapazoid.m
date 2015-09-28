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
nGrid = problem.options.trapazoid.nGrid;
if problem.options.verbose > 0
    fprintf('  -> Transcription via trapazoid method, nGrid = %d\n',nGrid);
end

%%%% Method-specific details to pass along to solver:

% Quadrature weights for trapazoid integration:
problem.func.weights = ones(nGrid,1);
problem.func.weights([1,end]) = 0.5;

% Trapazoid integration calculation of defects:
problem.func.defectCst = @computeDefects;

%%%% The key line - solve the problem by direct collocation:
soln = directCollocation(problem);

%%%% Method-specific interpolation
soln.interp.state = @(t)( interp1(tSoln',xSoln',t','linear',nan)' );
soln.interp.control = @(t)( interp1(tSoln',uSoln',t','linear',nan)' );

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function [defects, defectsGrad] = computeDefects(dt,x,f,dtGrad,xGrad,fGrad)
%
% TODO: Documentation
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

