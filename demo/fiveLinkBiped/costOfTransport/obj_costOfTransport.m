function [dObj, dObjGrad] = obj_costOfTransport(u,p)
% [dObj, dObjGrad] = obj_torqueSquared(u)
%
% This function computes the integrand of the cost of transport integral
%

empty = zeros(size(u(1,:)));

sn = u(6:10,:);
sp = u(11:15,:);

if nargout == 1 % numerical gradients
    
    dObj = autoGen_obj_costOfTransport(...
        p.m1,p.m2,p.m3,p.m4,p.m5,...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        sn(1,:),sn(2,:),sn(3,:),sn(4,:),sn(5,:),...
        sp(1,:),sp(2,:),sp(3,:),sp(4,:),sp(5,:),...
        p.g,p.stepLength,p.alpha,empty);
    
else  %Analytic gradients
    
    [dObj,fz,fzi] = autoGen_obj_costOfTransport(...
        p.m1,p.m2,p.m3,p.m4,p.m5,...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        sn(1,:),sn(2,:),sn(3,:),sn(4,:),sn(5,:),...
        sp(1,:),sp(2,:),sp(3,:),sp(4,:),sp(5,:),...
        p.g,p.stepLength,p.alpha,empty);
    
    dObjGrad = zeros(26,length(dObj));  % 26 = 1 + (5 + 5) + (5 + 10) = time + angle + rate + control + slack
    dObjGrad(fzi,:) = fz;
    
end

end