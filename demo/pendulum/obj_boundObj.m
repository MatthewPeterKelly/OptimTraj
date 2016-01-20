function [dObj, dObjGrad] = obj_boundObj(xF,xF_des)
% [dObj, dObjGrad] = obj_boundObj(u)
%
% This is a sample bound objective for testing analtyic gradients.

% weight on bound cost
Q = 100000*diag([1, 1]);

    % add cost for being away from final desired state.
    % NOTE: probably need to angle wrapping
    dObj = 0.5 * (xF-xF_des)' * Q * (xF-xF_des);
    
if nargout > 1  %Analytic gradients
  
    dObjGrad = zeros(1,6); %6 = [t0(1), x0(2), tF(1), xF(2)];
    
    dObjGrad(1,5:6) = (xF-xF_des)' * Q;  %gradient obj w.r.t. state xF
   
end

end