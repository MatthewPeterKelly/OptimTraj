function [c, ceq] = periodicGait(xF,x0,param)
% [c, ceq] = periodicGait(xF,x0,param)
%
% This function applies the non-linear boundary constraint to ensure that
% there gait is periodic, with the correct heel-strike map and left-right
% mapping.
%

% Unpack state
q = xF(1:5);
dq = xF(6:10);

% Compute heel-strike
[qNext, dqNext] = heelStrikeMap(q,dq,param);

% Pack up constraints
c = [];
ceq = x0 - [qNext; dqNext];

end