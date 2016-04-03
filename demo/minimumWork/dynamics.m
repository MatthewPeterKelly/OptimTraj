function [dx, dxGrad] = dynamics(x,u)
% [dx, dxGrad] = dynamics(x,u)
%
% Computes the dynamics (and gradients) for a 1d point mass on a
% friction-less plane with a force actuator.
%

% q = x(1,:);   %Position
dq = x(2,:);   %Velocity
ddq = u;  %Acceleration

dx = [dq;ddq];   %Pack up derivative of state

if nargout == 2   % Analytic gradients
    nTime = length(u);
    
    dqGrad = zeros(1,4,nTime); %4 = [time + pos + vel + force];
    dqGrad(1,3,:) = 1; %gradient dq wrt dq
    
    ddqGrad = zeros(1,4,nTime);  %4 = [time + angle + rate + torque];
    ddqGrad(1,4,:) = 1;  %gradient ddq wrt u
    
    dxGrad = cat(1, dqGrad, ddqGrad);
    
end

end