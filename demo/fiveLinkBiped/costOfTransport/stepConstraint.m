function [c, ceq, cGrad, ceqGrad] = stepConstraint(x0,xF,p)
% [c, ceq, cGrad, ceqGrad] = stepConstraint(x0,xF,p)
%
% This function applies the non-linear boundary constraint to ensure that
% there gait is periodic, with the correct heel-strike map and left-right
% mapping.
%

if nargout == 2 %Numerical gradients
    
    % Gait must be periodic
    ceq1 = cst_heelStrike(x0,xF,p);
    
    % Foot collision (toe-off and heel-strike) velocity
    c = cst_footVel(x0,xF,p);
    
    % Step length and height
    ceq2 = cst_stepLength(xF,p); 
    
    % Pack up equality constraints:
    ceq = [ceq1;ceq2];
       
    
else %Analytic gradients
    
    % Gait must be periodic
    [ceq1, ceqGrad1] = cst_heelStrike(x0,xF,p);
    
    % Foot collision (toe-off and heel-strike) velocity
    [c, cGrad] = cst_footVel(x0,xF,p);
    
    % Step length and height
    [ceq2, ceqGrad2] = cst_stepLength(xF,p);
    
        % Pack up equality constraints:
    ceq = [ceq1;ceq2];
    ceqGrad = [ceqGrad1;ceqGrad2];
    
end

end