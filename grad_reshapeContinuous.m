function grad = grad_reshapeContinuous(gradRaw,gradInfo)
% grad = grad_reshapeContinuous(gradRaw,gradInfo)
%
% TrajOpt utility function.
%
% This function converts the raw gradients from the user function into
% gradients with respect to the decision variables.
%
% INPUTS:
%   stateRaw = [nOutput,nInput,nTime]
%
% OUTPUTS:
%   grad = [nOutput,nTime,nDecVar]
%

if isempty(gradRaw)
    grad = [];
else
    [nOutput, ~, nTime] = size(gradRaw);
    
    grad = zeros(nOutput,nTime,gradInfo.nDecVar);
    
    % First, loop through and deal with time.
    timeGrad = gradRaw(:,1,:); timeGrad = permute(timeGrad,[1,3,2]);
    for iOutput=1:nOutput
        A = ([1;1]*timeGrad(iOutput,:)).*gradInfo.alpha;
        grad(iOutput,:,gradInfo.tIdx) = permute(A,[3,2,1]);
    end
    
    % Now deal with state and control:
    for iOutput=1:nOutput
        for iTime=1:nTime
            B = gradRaw(iOutput,2:end,iTime);
            grad(iOutput,iTime,gradInfo.xuIdx(:,iTime)) = permute(B,[3,1,2]);
        end
    end
end

end