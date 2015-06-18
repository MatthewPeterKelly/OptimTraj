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
%
%   nlpOpt = option struct to be passed to fmincon, created via optimset().
%
%   options = options for the transcription algorithm (this function)
%
%       .method = string to pick which method is used for transcription
%           'trapazoid'
%
%       .[method] = a struct to pass method-specific parameters. For
%       example, to pass the number of grid-points to the trapazoid method,
%       create a field .trapazoid.nGrid = [number of grid-points].
%
%
%       --TODO--
%       .verbose = integer
%           0 = no display
%           1 = default
%           2 = display warnings, overrides fmincon display setting
%           3 = debug
%
%
%
%
% OUTPUT: "soln"  --  struct with fields:
%
%

P = defaultTrajOpt(problem);

switch P.options.method
    case 'trapazoid'
        
        soln = trapazoid(P);
    otherwise
        error('Invalid method. Type: ''help trajOpt'' for a valid list.');
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
    
    % Check initial and final times:
    if ~checkField(problem.bounds,'initialTime')
        error('Field ''initialTime'' cannot be ommitted from ''problem.bounds''');
    end
    if ~checkField(problem.bounds,'finalTime'),
        error('Field ''finalTime'' cannot be ommitted from ''problem.bounds'''); end
    

    % Check to see if bounds exist, and fill in defaults if not
    if ~checkField(problem.bounds,'state')
        problem.bounds.state.low = -inf(nState,1);
        problem.bounds.state.upp = inf(nState,1);
    end 
    if ~checkField(problem.bounds,'initialState')
        problem.bounds.initialState.low = problem.bounds.state.low;
        problem.bounds.initialState.upp = problem.bounds.state.upp;
    end 
    if ~checkField(problem.bounds,'finalState')
        problem.bounds.finalState.low = problem.bounds.state.low;
        problem.bounds.finalState.upp = problem.bounds.state.upp;
    end 
    if ~checkField(problem.bounds,'control')
       problem.bounds.control.low = -inf(nControl,1);
       problem.bounds.control.upp = inf(nControl,1);
    end
    
    % Check the size (and existance) of .low and .upp
    checkLowUpp(problem.bounds.initialTime,1,1,'initialTime');
    checkLowUpp(problem.bounds.finalTime,1,1,'finalTime');
    checkLowUpp(problem.bounds.state,nState,1,'state');
    checkLowUpp(problem.bounds.initialState,nState,1,'initialState');
    checkLowUpp(problem.bounds.finalState,nState,1,'finalState');
    checkLowUpp(problem.bounds.control,nControl,1,'control');
    
end


%%%% Check options for trajOpt

if ~checkField(problem,'options'), problem.options = []; end

if ~checkField(problem.options, 'method')
    problem.options.method = 'trapazoid';   end

if ~checkField(problem.options, 'verbose')
    problem.options.verbose = 1;            end


%%%% Default options for fmincon
if ~checkField(problem,'nlpOpt')
   problem.nlpOpt = optimset('fmincon'); 
   problem.nlpOpt.Display = 'iter';
end

% override fmincon display for extreme verbose options
if problem.options.verbose == 0
    problem.nlpOpt.Display = 'off';
elseif problem.options.verbose == 4
    problem.nlpOpt.Display = 'iter-detailed';
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


function checkLowUpp(input,nRow,nCol,name)
%
% This function checks that input has the following is true:
%   size(input.low) == [nRow, nCol]
%   size(input.upp) == [nRow, nCol]

if ~checkField(input,'low')
    error(['Field ''low'' cannot be ommitted from problem.bounds.' name '']);
end

if ~checkField(input,'upp')
    error(['Field ''lupp'' cannot be ommitted from problem.bounds.' name '']);
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

end