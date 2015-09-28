function soln = gpopsWrapper(problem)
% soln = gpopsWrapper(problem)
%
% This function is a wrapper that converts the standard input for trajOpt
% into a call to GPOPS2, a commercially available transcription software
% for matlab. You can purchase and download it at http://www.gpops2.com/
%
% GPOPS2 implements an adaptive transcription method - it adjusts both the
% number of trajectory segments and the order of the interpolating
% polynomial in each segment. Most GPOPS features are available through
% this interface, with the exception of multiple phase problems.
%

% Print out some solver info if desired:
if problem.options.verbose > 0
   disp('Transcription using GPOPS2');
end





disp('here');

%%%% Method-specific interpolation
soln.interp.state = @(t)( interp1(tSoln',xSoln',t','linear',nan)' );
soln.interp.control = @(t)( interp1(tSoln',uSoln',t','linear',nan)' );

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function [defects, defectsGrad] = computeDefects(dt,x,f,dtGrad,xGrad,fGrad)
%
% This function computes the defects that are used to enforce the
% continuous dynamics of the system along the trajectory.
%
% INPUTS:
%   dt = time step (scalar)
%   x = [nState, nTime] = state at each grid-point along the trajectory
%   f = [nState, nTime] = dynamics of the state along the trajectory
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   dtGrad = [2,1] = gradient of time step with respect to [t0; tF]
%   xGrad = [nState,nTime,nDecVar] = gradient of trajectory wrt dec vars
%   fGrad = [nState,nTime,nDecVar] = gradient of dynamics wrt dec vars
%
% OUTPUTS:
%   defects = [nState, nTime-1] = error in dynamics along the trajectory
%   defectsGrad = [nState, nTime-1, nDecVars] = gradient of defects
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

