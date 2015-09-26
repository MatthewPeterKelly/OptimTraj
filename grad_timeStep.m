function [dt, dtGrad] = grad_timeStep(t,gradInfo)
%
% TrajOpt utility function
%
% Computes the time step and its gradient
%
% dt = [1,1]
% dtGrad = [1,nz]
%

nTime = length(t);

dt = (t(end)-t(1))/(nTime-1);
dtGrad = zeros(1,gradInfo.nDecVar);

dtGrad(1,gradInfo.tIdx(1)) = -1/(nTime-1);
dtGrad(1,gradInfo.tIdx(2)) = 1/(nTime-1);

end