function [obj, objGrad] = pathObjectiveTest(t,x,u)
% [obj, objGrad] = pathObjectiveTest(t,x,u)
%
% Computes a test objective function with time, state, and control
% dependence.
% 


obj = u.^2 + x(2,:).^2 + t;

if nargout == 2  % Analytic gradients
    nTime = length(u);
    
    objGrad = zeros(4,nTime); %4 = [time + angle + rate + torque];
    
    objGrad(1,:) = ones(size(t));
    objGrad(3,:) = 2*x(2,:);
    objGrad(4,:) = 2*u;  %gradient obj wrt u
    
end

end