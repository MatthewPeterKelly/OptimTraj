function [obj, objGrad] = pathObjective_dev(t,x,u)
%  [obj, objGrad] = pathObjective_dev(t,x,u)
%
% Computes the objective function (and gradients) for the simple pendulum
% with a complicated cost function that has time, state, and control
% dependency. This is primarily used for checking gradients.
%

obj = u.^2 + x(2,:).^2 + t;

if nargout == 2  % Analytic gradients
    nTime = length(u);
    
    objGrad = zeros(4,nTime); %4 = [time + angle + rate + torque];
    
    objGrad(1,:) = 1;
    objGrad(3,:) = 2*x(2,:);
    objGrad(4,:) = 2*u;  %gradient obj wrt u
    
end

end