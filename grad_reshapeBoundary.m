function CC = grad_reshapeBoundary(C,gradInfo)
%
% This function takes a boundary constraint or objective from the user
% and expands it to match the full set of decision variables
%

CC = zeros(size(C,1),gradInfo.nDecVar);
CC(:,gradInfo.bndIdxMap) = C;

end