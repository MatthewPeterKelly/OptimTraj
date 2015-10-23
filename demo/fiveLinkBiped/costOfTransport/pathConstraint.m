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


if nargout == 2 % numerical gradients
    
    c = [...
        q1-q2;    %Stance knee joint limit
        q5-q4];   %Swing knee joint limit
    
else %Analytic gradients
    
    
    %%%% Joint Limits
    c = [...
        q1-q2;    %Stance knee joint limit
        q5-q4];   %Swing knee joint limit
       
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
    
    %%%% Slack Variables
    dq = x(6:10<
    [ceq, ceqZ, ceqZi] = autoGen_cst_costOfTransport(...
        dq1,dq2,dq3,dq4,dq5,...
        u1,u2,u3,u4,u5,....
        sn1,sn2,sn3,sn4,sn5,...
        sp1,sp2,sp3,sp4,sp5,...
        empty);
    
    
    ceqGrad = [];
end

end