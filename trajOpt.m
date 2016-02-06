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
%       Analytic Gradients:
%           Both the "trapezoid" and "hermiteSimpson" methods in TrajOpt
%       support analytic gradients. Type help "trapezoid" or help
%       "hermiteSimpson" for details on how to pass gradients to the
%       transcription method.
%
%
%   bounds - struct with bounds for the problem:
%
%       .initialTime.low = [scalar]
%       .initialTime.upp = [scalar]
%
%       .finalTime.low = [scalar]
%       .finalTime.upp = [scalar]
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
%           'trapezoid'
%           'hermiteSimpson'
%           'chebyshev'
%           'rungeKutta'
%           'gpops'     ( Must have GPOPS2 installed )
%
%       .[method] = a struct to pass method-specific parameters. For
%       example, to pass the number of grid-points to the trapezoid method,
%       create a field .trapezoid.nGrid = [number of grid-points].
%
%       .verbose = integer
%           0 = no display
%           1 = default
%           2 = display warnings, overrides fmincon display setting
%           3 = debug
%
%       .defaultAccuracy = {'low','medium','high'}
%           Sets the default options for each transcription method
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
%   * Both the trapezoid and hermiteSimpson method support additional
%   features. The first is that they can use analytic gradients, if
%   provided by the user. ">> help trapezoid" or ">> help hermiteSimpson"
%   for more information. These methods additionally provide error
%   estimates in two forms. The continuous collocation constraint error is
%   provided as an additional function handle (soln.interp.collCst), and
%   the absolute local error estimate for each segment is provided in
%   soln.info.error.
%

problem = inputValidation(problem);   %Check inputs
problem = getDefaultOptions(problem); % Complete options struct

% Loop over the options struct to solve the problem
nIter = length(problem.options);
soln(nIter) = struct('grid',[],'interp',[],'info',[],'problem',[]); 
P = problem;  %Temp variable for passing on each iteration
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
        case 'trapezoid'
            soln(iter) = trapezoid(P);
        case 'hermiteSimpson'
            soln(iter) = hermiteSimpson(P);
        case 'chebyshev'
            soln(iter) = chebyshev(P);
        case 'multiCheb'
            soln(iter) = multiCheb(P);
        case 'rungeKutta'
            soln(iter) = rungeKutta(P);
        case 'gpops'
            soln(iter) = gpopsWrapper(P);
        otherwise
            error('Invalid method. Type: ''help trajOpt'' for a valid list.');
    end
    
end
end