function [c, ceq, cGrad, ceqGrad] = pathConstraint(x)
% [c, ceq, cGrad, ceqGrad] = pathConstraint(x)
%
% This function implements a simple path constraint to keep the knee joint
% of the robot from hyer-extending.
%

q1 = x(1,:);
q2 = x(2,:);
q4 = x(4,:);
q5 = x(5,:);

c = [...
    q1-q2;    %Stance knee joint limit
    q5-q4];   %Swing knee joint limit

ceq = [];

if nargout == 4 %Analytic gradients
    % Gradients with respect to:
    % [t,q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5,u1,u2,u3,u4,u5] = 1+5+5+5
    nCst = 2;   %stance leg ; swing leg
    nGrad = 16;  %time, angles, rates, torques
    nTime = size(x,2);
    cGrad = zeros(nCst,nGrad,nTime);
    cGrad(1,3,:) = -1;  % cst stance wrt q2
    cGrad(1,2,:) = 1; % cst stance wrt q1
    cGrad(2,5,:) = -1;  % cst swing wrt q4
    cGrad(2,6,:) = 1; % cst swing wrt q5
    
    ceqGrad = [];
end

end