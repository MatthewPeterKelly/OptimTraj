function [obj, objGrad] = pathObjective(u)
% [obj, objGrad] = pathObjective(u)
%
% Computes the objective function (and gradients) for the simple pendulum
%

obj = u.^2;

if nargout == 2  % Analytic gradients
    nTime = length(u);
    
    objGrad = zeros(4,nTime); %4 = [time + angle + rate + torque];
    
    objGrad(4,:) = 2*u;  %gradient obj wrt u
    
end

end