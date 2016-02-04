function [obj, objGrad] = DEV_pathObjective(t,x,u)
%  [obj, objGrad] = DEV_pathObjective(t,x,u)
%
% Computes the objective function (and gradients) for the simple pendulum
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