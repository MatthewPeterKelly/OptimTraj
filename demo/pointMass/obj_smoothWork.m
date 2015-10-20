function [obj, objGrad] = obj_smoothWork(x,u,alpha, beta)
% [obj, objGrad] = obj_smoothWork(x,u,alpha,beta)
%
% Computes the objective function (and gradients) for the simple pendulum
%
% x = state = [pos;vel]
% u = force
% alpha = smoothing parameter for abs()
% beta = addative torque-squared smoothing

% q = x(1,:);   %Position
dq = x(2,:);   %Velocity
power = dq.*u;  % power = speed * force

% work = smoothing + integral(abs(power)) 
obj = beta*u.^2 + power.*tanh(power/alpha);  %exponential smoothing of abs(power)

if nargout == 2  % Analytic gradients
    nTime = length(u);
    
    objGrad = zeros(4,nTime); %4 = [time + angle + rate + torque];
    
    % dObj/dDq
    objGrad(3,:) = u.*tanh((dq.*u)./alpha) - (dq.*u.^2.*(tanh((dq.*u)./alpha).^2 - 1))./alpha;
    
    % dObj/dU
    objGrad(4,:) = 2*beta*u + dq.*tanh((dq.*u)./alpha) - (dq.^2.*u.*(tanh((dq.*u)./alpha).^2 - 1))./alpha;  
    
end

end