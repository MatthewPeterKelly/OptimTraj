function [c, ceq, cGrad, ceqGrad] = cstSlackPower(x,u)
% [c, ceq, cGrad, ceqGrad] = cstSlackPower(x,u)
%
% This is the path constraint that is used to ensure that the sum of the
% difference between the two slack variables is the power.
%
% See Bett's Book, chapter 1.16 for mathematical details.
%

% q = x(1,:);   %Position
dq = x(2,:);   %Velocity
f = u(1,:);  %Force
power = dq.*f;  % power = speed * force

% Slack variables
s1 = u(2,:);
s2 = u(3,:);

c = [];
ceq = power - (s1-s2);

if nargout > 2   %Analytic gradients
    
    nTime = length(f);
    
    cGrad = [];
    
    ceqGrad = zeros(1,6,nTime); %6 = [time + angle + rate + torque + slack1 + slack 2];
    
    ceqGrad(1,3,:) = f; % derivative wrt rate
    ceqGrad(1,4,:) = dq; % derivative wrt force
    ceqGrad(1,5,:) = -1; % derivative wrt s1
    ceqGrad(1,6,:) = 1; % derivative wrt s2
    
end

end