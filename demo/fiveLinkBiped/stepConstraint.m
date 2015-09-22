function [c, ceq] = stepConstraint(xF,x0,param)
% [c, ceq] = stepConstraint(xF,x0,param)
%
% This function applies the non-linear boundary constraint to ensure that
% there gait is periodic, with the correct heel-strike map and left-right
% mapping.
%

stepLength = 0.4;   %Desired step length in meters
slope = 0.0;

% Unpack state
q = xF(1:5);
dq = xF(6:10);

% Compute heel-strike
[qNext, dqNext] = heelStrikeMap(q,dq,param);
defects = x0 - [qNext; dqNext];

% Compute the distance moved by the center of mass:
G = getComPos([qNext,q],param);
delGx = G(1,2) - G(1,1);  %Change in x position
delGy = G(2,2) - G(2,1);  %Change in y position
stepHeight = slope/stepLength;
stepCst = [delGx - stepLength; delGy - stepHeight];

% Pack up constraints
c = [];
ceq = [defects; stepCst];

end