function soln = trajOpt(problem)
% soln = trajOpt(problem)
%
% Solves a trajectory optimization problem.
%
% INPUT: "problem" -- struct with fields:
%
%   func -- struct for user-defined functions, passed as function handles
%
%       Input Notes:
%               t = [1, nTime] = time vector (grid points)
%               x = [nState, nTime] = state vector at each grid point
%               u = [nControl, nTime] = control vector at each grid point
%               t0 = scalar = initial time
%               tF = scalar = final time
%               x0 = [nState, 1] = initial state
%               xF = [nState, 1] = final state
%
%       dx = dynamics(t,x,u)
%               dx = [nState, nTime] = dx/dt = derivative of state wrt time
%
%       dObj = pathObj(t,x,u)
%               dObj = [1, nTime] = integrand from the cost function
%
%       obj = bndObj(t0,x0,tF,xF)
%               obj = scalar = objective function for boundry points
%
%       [c, ceq] = pathCst(t,x,u)
%               c = column vector of inequality constraints  ( c <= 0 )
%               ceq = column vector of equality constraints ( c == 0 )
%
%       [c, ceq] = bndCst(t0,x0,tF,xF)
%               c = column vector of inequality constraints  ( c <= 0 )
%               ceq = column vector of equality constraints ( c == 0 )
%
%       How to pass parameters to your functions:
%           - suppose that your dynamics function is pendulum.m and it
%           accepts a struct of parameters p. When you are setting up the
%           problem, define the struc p in your workspace and then use the
%           following command to pass the function:
%               problem.func.dynamics = @(t,x,u)( pendulum(t,x,u,p) );
%
%   bounds - struct with bounds for the problem:
%
%       initialTime.low = [1, 1]
%       initialTime.upp = [1, 1]
%
%       finalTime.low = [1, 1]
%       finalTime.upp = [1, 1]
%
%       .state.low = [nState,1] = lower bound on the state
%       .state.upp = [nState,1] = lower bound on the state
%
%       .initialState.low = [nState,1]
%       .initialState.upp = [nState,1]
%
%       .finalState.low = [nState,1]
%       .finalState.upp = [nState,1]
%
%       .control.low = [nControl, 1]
%       .control.upp = [nControl, 1]
%
%
%
%   guess - struct with an initial guess at the trajectory
%
%       .time = [1, nGridGuess]
%       .state = [nState, nGridGuess]
%       .control = [nControl, nGridGuess]
%
%   options = options for the transcription algorithm (this function)
%
%       .nlpOpt = options to pass through to fmincon
%
%       .method = string to pick which method is used for transcription
%           'trapazoid'
%
%       .[method] = a struct to pass method-specific parameters. For
%       example, to pass the number of grid-points to the trapazoid method,
%       create a field .trapazoid.nGrid = [number of grid-points].
%
%       .verbose = integer
%           0 = no display
%           1 = default
%           2 = display warnings, overrides fmincon display setting
%           3 = debug
%
%       * if options is a struct array, the trajOpt will run the optimization
%       by running options(1) and then using the result to initialize a new
%       solve with options(2) and so on, until it runs options (end). This
%       allows for successive grid and tolerance opdates.
%
%
%
%
%
% OUTPUT: "soln"  --  struct with fields:
%
%   .grid = trajectory at the grid-points used by the transcription method
%       .time = [1, nTime]
%       .state = [nState, nTime]
%       .control = [nControl, nTime];
%
%   .interp = functions for interpolating state and control for arbitrary
%       times long the trajectory. The interpolation method will match the
%       underlying transcription method. This is particularily important
%       for high-order methods, where linear interpolation between the
%       transcription grid-points will lead to large errors. If the
%       requested time is not on the trajectory, the interpolation will
%       return NaN.
%
%       .state = @(t) = given time, return state
%           In: t = [1,n] vector of time
%           Out: x = [nState,n] state vector at each point in time
%
%       .control = @(t) = given time, return control
%           In: t = [1,n] vector of time
%           Out: u = [nControl,n] state vector at each point in time
%
%   .info = information about the optimization run
%       .nlpTime = time (seconds) spent in fmincon
%       .exitFlag = fmincon exit flag
%       .objVal = value of the objective function
%       .[all fields in the fmincon "output" struct]
%
%   .problem = the problem as it was passed to the low-level transcription,
%       including the all default values that were used
%
%
%   * If problem.options was a struct array, then soln will also be a
%   struct array, with soln(1) being the solution on the first iteration,
%   and soln(end) being the final solution.
%

problem = defaultTrajOpt(problem);
P = problem; P.options = [];

% Loop over the options struct to solve the problem
nIter = length(problem.options);
soln(nIter) = struct('grid',[],'interp',[],'info',[],'problem',[]);  %Initialize struct array
for iter=1:nIter
    P.options = problem.options(iter);
    
    if P.options.verbose > 0    %then print out iteration count:
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
        disp(['Running TrajOpt, iteration ' num2str(iter)]);
    end
    
    if iter > 1  %Use previous soln as new guess
        P.guess = soln(iter-1).grid;
    end
    
    %%%% This is the key part: call the underlying transcription method:
    switch P.options.method
        case 'trapazoid'
            soln(iter) = trapazoid(P);
        case 'chebyshev'
            soln(iter) = chebyshev(P);
        otherwise
            error('Invalid method. Type: ''help trajOpt'' for a valid list.');
    end
    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       SUB-FUNCTIONS                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function problem = defaultTrajOpt(problem)
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

if ~checkField(problem,'func')
    error('Field ''func'' cannot be ommitted from ''problem''');
else
    if ~checkField(problem.func,'dynamics')
        error('Field ''dynamics'' cannot be ommitted from ''problem.func'''); end
    if ~checkField(problem.func,'pathObj'), problem.func.pathObj = []; end
    if ~checkField(problem.func,'bndObj'), problem.func.bndObj = []; end
    if ~checkField(problem.func,'pathCst'), problem.func.pathCst = []; end
    if ~checkField(problem.func,'bndCst'), problem.func.bndCst = []; end
end

%%%% Check the initial guess (also compute nState and nControl):
if ~checkField(problem, 'guess')
    error('Field ''guess'' cannot be ommitted from ''problem''');
else
    if ~checkField(problem.guess,'time')
        error('Field ''time'' cannot be ommitted from ''problem.guess'''); end
    if ~checkField(problem.guess, 'state')
        error('Field ''state'' cannot be ommitted from ''problem.guess'''); end
    if ~checkField(problem.guess, 'control')
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
if ~checkField(problem,'bounds')
    error('Field ''bounds'' cannot be ommitted from ''problem''');
else
    
    if ~checkField(problem.bounds,'initialTime')
        problem.bounds.initialTime = []; end
    problem.bounds.initialTime = ...
        checkLowUpp(problem.bounds.initialTime,1,1,'initialTime');
    
    if ~checkField(problem.bounds,'finalTime')
        problem.bounds.finalTime = []; end
    problem.bounds.finalTime = ...
        checkLowUpp(problem.bounds.finalTime,1,1,'finalTime');
    
    if ~checkField(problem.bounds,'state')
        problem.bounds.state = []; end
    problem.bounds.state = ...
        checkLowUpp(problem.bounds.state,nState,1,'state');
    
    if ~checkField(problem.bounds,'initialState')
        problem.bounds.initialState = []; end
    problem.bounds.initialState = ...
        checkLowUpp(problem.bounds.initialState,nState,1,'initialState');
    
    if ~checkField(problem.bounds,'finalState')
        problem.bounds.finalState = []; end
    problem.bounds.finalState = ...
        checkLowUpp(problem.bounds.finalState,nState,1,'finalState');
    
    if ~checkField(problem.bounds,'control')
        problem.bounds.control = []; end
    problem.bounds.control = ...
        checkLowUpp(problem.bounds.control,nControl,1,'control');
    
end


%%%% Check options for trajOpt

if ~checkField(problem,'options')
    problem.options = struct('method',[],'verbose',[],'nlpOpt',[]);
end

% Loop over each iteration:
for iter = 1:length(problem.options)
    
    if ~checkField(problem.options(iter), 'method')
        problem.options(iter).method = 'trapazoid';   end
    
    if ~checkField(problem.options(iter), 'verbose')
        problem.options(iter).verbose = 1;            end
    
    
    %%%% Default options for fmincon
    if ~checkField(problem.options(iter),'nlpOpt')
        problem.options(iter).nlpOpt = optimset('fmincon');
        problem.options(iter).nlpOpt.Display = 'iter';
    end
    
    % override fmincon display for extreme verbose options
    if problem.options(iter).verbose == 0
        problem.options(iter).nlpOpt.Display = 'off';
    elseif problem.options(iter).verbose == 4
        problem.options(iter).nlpOpt.Display = 'iter-detailed';
    end
    
end

end


function valid = checkField(input,field)
%
% This function returns true if the field exists and is non-empty
%

valid = false;   %Assume the worst
if isfield(input,field)   %Check that field exists
    if ~isempty(input.(field))   %Check that it's non-empty
        valid = true;  %yay! A valid field
    end
end

end


function input = checkLowUpp(input,nRow,nCol,name)
%
% This function checks that input has the following is true:
%   size(input.low) == [nRow, nCol]
%   size(input.upp) == [nRow, nCol]

if ~checkField(input,'low')
    input.low = -inf(nRow,nCol);
end

if ~checkField(input,'upp')
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