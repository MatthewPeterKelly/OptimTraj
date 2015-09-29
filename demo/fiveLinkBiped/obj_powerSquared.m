function [dObj, dObjGrad] = obj_powerSquared(x,u)
% [dObj, dObjGrad] = obj_powerSquared(x,u)
%
% This function computes the integrand of the cost of transport integral
%

dq = x(6:10,:);

if nargout == 1 % numerical gradients
    
    dObj = autoGen_obj_powerSquared(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:));
    
else  %Analytic gradients
    
    [dObj,fz,fzi] = autoGen_obj_powerSquared(...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:));
    dObjGrad = zeros(16,length(dObj));  % 16 = 1 + 5 + 5 + 5 = time + angle + rate + control
    dObjGrad(fzi,:) = fz;
    
end

end