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
%
%       --TODO-- 
%  
%
%
%
% OUTPUT: "soln"  --  struct with fields:
%
%   

P = problem;

% --TODO--
%   Fill in default values where necessary
%

switch P.options.method
    case 'trapazoid'
        soln = trapazoid(problem);
    otherwise
        error('Invalid method. Type: ''help trajOpt'' for a valid list.');
end

end