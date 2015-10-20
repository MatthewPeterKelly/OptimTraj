function [dx, dxGrad] = cstDyn(x,u)
% [dx, dxGrad] = cstDyn(x,u)
%
% Computes the dynamics (and gradients) for a 1d point mass on a
% friction-less plane with a force actuator.
%

% q = x(1,:);   %Position
dq = x(2,:);   %Velocity
ddq = u(1,:);  %Acceleration

dx = [dq;ddq];   %Pack up derivative of state

if nargout == 2   % Analytic gradients
    nTime = length(u);
    
    dqGrad = zeros(1,6,nTime); %6 = [time + pos + vel + force + slack];
    dqGrad(1,3,:) = 1; %gradient dq wrt dq
    
    ddqGrad = zeros(1,6,nTime);  %6 = [time + angle + rate +  force + slack];
    ddqGrad(1,4,:) = 1;  %gradient ddq wrt u
    
    dxGrad = cat(1, dqGrad, ddqGrad);
    
end

end