function problem = inputValidation(problem)
%
% This function runs through the problem struct and sets any missing fields
% to the default value. If a mandatory field is missing, then it throws an
% error.
%
% INPUTS:
%   problem = a partially completed problem struct
%
% OUTPUTS:
%   problem = a complete problem struct, with validated fields
%


%%%% Check the function handles:

if ~isfield(problem,'func')
    error('Field ''func'' cannot be ommitted from ''problem''');
else
    if ~isfield(problem.func,'dynamics')
        error('Field ''dynamics'' cannot be ommitted from ''problem.func'''); end
    if ~isfield(problem.func,'pathObj'), problem.func.pathObj = []; end
    if ~isfield(problem.func,'bndObj'), problem.func.bndObj = []; end
    if ~isfield(problem.func,'pathCst'), problem.func.pathCst = []; end
    if ~isfield(problem.func,'bndCst'), problem.func.bndCst = []; end
end

%%%% Check the initial guess (also compute nState and nControl):
if ~isfield(problem, 'guess')
    error('Field ''guess'' cannot be ommitted from ''problem''');
else
    if ~isfield(problem.guess,'time')
        error('Field ''time'' cannot be ommitted from ''problem.guess'''); end
    if ~isfield(problem.guess, 'state')
        error('Field ''state'' cannot be ommitted from ''problem.guess'''); end
    if ~isfield(problem.guess, 'control')
        error('Field ''control'' cannot be ommitted from ''problem.guess'''); end
    
    % Compute the size of the time, state, and control based on guess
    [checkOne, nTime] = size(problem.guess.time);
    [nState, checkTimeState] = size(problem.guess.state);
    [nControl, checkTimeControl] = size(problem.guess.control);
    
    if nTime < 2 || checkOne ~= 1
        error('guess.time must have dimensions of [1, nTime], where nTime > 1');
    end
    
    if checkTimeState ~= nTime
        error('guess.state must have dimensions of [nState, nTime]');
    end
    if checkTimeControl ~= nTime
        error('guess.control must have dimensions of [nControl, nTime]');
    end
    
end

%%%% Check the problem bounds:
if ~isfield(problem,'bounds')
    error('Field ''bounds'' cannot be ommitted from ''problem''');
else
    
    if ~isfield(problem.bounds,'initialTime')
        problem.bounds.initialTime = []; end
    problem.bounds.initialTime = ...
        checkLowUpp(problem.bounds.initialTime,1,1,'initialTime');
    
    if ~isfield(problem.bounds,'finalTime')
        problem.bounds.finalTime = []; end
    problem.bounds.finalTime = ...
        checkLowUpp(problem.bounds.finalTime,1,1,'finalTime');
    
    if ~isfield(problem.bounds,'state')
        problem.bounds.state = []; end
    problem.bounds.state = ...
        checkLowUpp(problem.bounds.state,nState,1,'state');
    
    if ~isfield(problem.bounds,'initialState')
        problem.bounds.initialState = []; end
    problem.bounds.initialState = ...
        checkLowUpp(problem.bounds.initialState,nState,1,'initialState');
    
    if ~isfield(problem.bounds,'finalState')
        problem.bounds.finalState = []; end
    problem.bounds.finalState = ...
        checkLowUpp(problem.bounds.finalState,nState,1,'finalState');
    
    if ~isfield(problem.bounds,'control')
        problem.bounds.control = []; end
    problem.bounds.control = ...
        checkLowUpp(problem.bounds.control,nControl,1,'control');
    
end

end


function input = checkLowUpp(input,nRow,nCol,name)
%
% This function checks that input has the following is true:
%   size(input.low) == [nRow, nCol]
%   size(input.upp) == [nRow, nCol]

if ~isfield(input,'low')
    input.low = -inf(nRow,nCol);
end

if ~isfield(input,'upp')
    input.upp = inf(nRow,nCol);
end

[lowRow, lowCol] = size(input.low);
if lowRow ~= nRow || lowCol ~= nCol
    error(['problem.bounds.' name ...
        '.low must have size = [' num2str(nRow) ', ' num2str(nCol) ']']);
end

[uppRow, uppCol] = size(input.upp);
if uppRow ~= nRow || uppCol ~= nCol
    error(['problem.bounds.' name ...
        '.upp must have size = [' num2str(nRow) ', ' num2str(nCol) ']']);
end

if sum(sum(input.upp-input.low < 0))
    error(...
        ['problem.bounds.' name '.upp must be >= problem.bounds.' name '.low!']);
end

end