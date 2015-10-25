function [c, ceq, cGrad, ceqGrad] = pathConstraint(state,control,param)
% [c, ceq, cGrad, ceqGrad] = pathConstraint(state,control,param)
%
% This function implements a simple path constraint to keep the knee joint
% of the robot from hyper-extending.
%

q1 = state(1,:);
q2 = state(2,:);
q4 = state(4,:);
q5 = state(5,:);

dq = state(6:10,:);
u = state(11:15,:);
sn = control(6:10,:);
sp = control(11:15,:);
empty = zeros(size(q1));

if nargout == 2 % numerical gradients
    
    y = autoGen_cst_swingFootHeight(...
        q1,q2,q4,q5,...
        param.l1,param.l2,param.l4,param.l5,param.stepLength,param.stepHeight);
    
    c = [...
        q1-q2;    %Stance knee joint limit
        q5-q4;   %Swing knee joint limit
        y];   %Swing foot height
        
     ceq = autoGen_cst_costOfTransport(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        sn(1,:),sn(2,:),sn(3,:),sn(4,:),sn(5,:),...
        sp(1,:),sp(2,:),sp(3,:),sp(4,:),sp(5,:),...
        empty);
    
else %Analytic gradients
    
    %%%% Swing foot clearance:
     [y,yz,yzi] = autoGen_cst_swingFootHeight(...
        q1,q2,q4,q5,...
        param.l1,param.l2,param.l4,param.l5,param.stepLength,param.stepHeight);
    
    %%%% Joint Limits
    c = [...
        q1-q2;    %Stance knee joint limit
        q5-q4;   %Swing knee joint limit
        y];   % Swing foot height
        
    % Gradients with respect to:
    % [t,q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5,u1,u2,u3,u4,u5] = 1+5+5+5
    nCst = 3;   %stance leg ; swing leg
    nGrad = 31;  %time, angles, rates, torques, torqueRate, slack
    nTime = size(state,2);
    cGrad = zeros(nCst,nGrad,nTime);
    cGrad(1,3,:) = -1;  % cst stance wrt q2
    cGrad(1,2,:) = 1; % cst stance wrt q1
    cGrad(2,5,:) = -1;  % cst swing wrt q4
    cGrad(2,6,:) = 1; % cst swing wrt q5
        
    %%%% Slack Variables
    [ceq, ceqZ, ceqZi] = autoGen_cst_costOfTransport(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        sn(1,:),sn(2,:),sn(3,:),sn(4,:),sn(5,:),...
        sp(1,:),sp(2,:),sp(3,:),sp(4,:),sp(5,:),...
        empty);
    
    ceqGrad = zeros(5,nGrad,nTime);  % 5 = number of slack constraints
    tmpHeight = zeros(1,1,nGrad);   %foot height
    tmpSlack = zeros(5,1,nGrad);  % Slack variables
    for i=1:nTime
        tmpHeight(yzi) = yz(:,i);
        cGrad(3,:,i) = tmpHeight(:,1,:);
        
        tmpSlack(ceqZi) = ceqZ(:,i);
        ceqGrad(:,:,i) = tmpSlack(:,1,:);         
    end
     
end

end