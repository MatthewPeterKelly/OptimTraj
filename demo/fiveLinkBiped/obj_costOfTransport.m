function [dObj, dObjGrad] = obj_costOfTransport(x,u,p)
% [dObj, dObjGrad] = obj_torqueSquared(u)
%
% This function computes the integrand of the cost of transport integral
%

dq = x(6:10,:);

if nargout == 1 % numerical gradients
    
    dObj = autoGen_obj_costOfTransport(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        p.m1, p.m2, p.m3, p.m4, p.m5, p.g,...
        p.alpha,p.stepLength);
    
else  %Analytic gradients
    
    [dObj,fz,fzi] = autoGen_obj_costOfTransport(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        p.m1, p.m2, p.m3, p.m4, p.m5, p.g,...
        p.alpha,p.stepLength);
    dObjGrad = zeros(16,length(dObj));  % 16 = 1 + 5 + 5 + 5 = time + angle + rate + control
    dObjGrad(fzi,:) = fz;
    
end

end