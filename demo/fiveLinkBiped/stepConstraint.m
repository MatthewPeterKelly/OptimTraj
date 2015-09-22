function [c, ceq] = stepConstraint(t0,x0,tF,xF,config)
% [c, ceq] = stepConstraint(t0,x0,tF,xF,config)
%
% This function applies the non-linear boundary constraint to ensure that
% there gait is periodic, with the correct heel-strike map and left-right
% mapping.
%

if config.slope ~= 0
    error('Non-zero walking slopes are not yet implemented!');
end


% Unpack state
q = xF(1:5);
dq = xF(6:10);

% Compute heel-strike
[qNext, dqNext] = heelStrikeMap(q,dq,config.param);
defects = x0 - [qNext; dqNext];

%%%% Compute non-linear constraints on the step:
[P5, dP5] = swingFootKinematics([qNext,q],[dqNext,dq],config.param);
% p0x = P5(1,1);  %Initial horizontal position of swing foot
% p0y = P5(2,1);  %Initial vertical position of swing foot
pFx = P5(1,2);  %Final horizontal position of swing foot
pFy = P5(2,2);  %Final vertical position of swing foot
% dp0x = dP5(1,1);  %Initial horizontal position of swing foot
dp0y = dP5(2,1);  %Initial vertical position of swing foot
% dpFx = dP5(1,2);  %Final horizontal position of swing foot
dpFy = dP5(2,2);  %Final vertical position of swing foot

% Step Length:
cstStepLength = pFx - config.stepLength;

% Step Height
cstStepHeight = pFy - 0.0;

% Foot leaves ground at start of step:
cstToeOff = -dp0y;

% Foot strikes ground at end of step:
cstHeelStrike = dpFy;

% Pack up constraints
c = [cstToeOff; cstHeelStrike];
ceq = [defects; cstStepLength; cstStepHeight];

end