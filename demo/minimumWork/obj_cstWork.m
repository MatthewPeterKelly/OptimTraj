function [obj, objGrad] = obj_cstWork(u)
% [obj, objGrad] = obj_smoothWork(u)
%
% Computes the objective function (and gradients) for the point mass
%

% Slack variables
s1 = u(2,:);
s2 = u(3,:);

% minimize the abs(power), where power = s1-s1, and s1 > 0, s2 > 0
obj = s1 + s2;

if nargout == 2  % Analytic gradients
    nTime = length(u);
    
    objGrad = zeros(6,nTime); %6 = [time + angle + rate + torque + slack1 + slack 2];
    
    % dObj/dS1
    objGrad(5,:) = 1;
    
    % dObj/dS2
    objGrad(6,:) = 1;  
    
end

end