function [dObj, dObjGrad] = boundObjective(xF,xF_target)
% [dObj, dObjGrad] = boundObjective(xF,xF_target)
%
% This is a sample bound objective for testing analtyic gradients.

% weight on bound cost
Q = diag([1, 1]);

    % add cost for being away from final desired state.
    % NOTE: probably need to angle wrapping
    dObj = 0.5 * (xF-xF_target)' * Q * (xF-xF_target);
    
if nargout > 1  %Analytic gradients
  
    dObjGrad = zeros(1,6); %6 = [t0(1), x0(2), tF(1), xF(2)];
    
    dObjGrad(1,5:6) = (xF-xF_target)' * Q;  %gradient obj w.r.t. state xF
   
end

end