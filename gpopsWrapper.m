function soln = gpopsWrapper(problem)
% soln = gpopsWrapper(problem)
%
% This function is a wrapper that converts the standard input for trajOpt
% into a call to GPOPS2, a commercially available transcription software
% for matlab. You can purchase and download it at http://www.gpops2.com/
%
% GPOPS2 implements an adaptive transcription method - it adjusts both the
% number of trajectory segments and the order of the interpolating
% polynomial in each segment. Many GPOPS features are available in TrajOpt,
% but not all. Notably, TrajOpt cannot solve multi-phase problems.
%
% Set any special GPOPS options by storing the 'setup' sturuct in the 
% problem.options.gpops struct.
%
% If using SNOPT, be careful about any constant terms in your constraints.
% When using numerical gradients, SNOPT drops any constant terms in your
% constraints, which is why it has non-zero bounds. This is exactly the
% opposite of the convention that FMINCON uses, where all constraint bounds
% must be zero. If your constraints have non-zero bounds, and you would
% like to use GPOPS with SNOPT as the solver, then manually set the fields
% in problem.gpops.phase.bounds.path and problem.gpops.eventgroup to
% include these bounds, and then remove them from the constraint function.
%


% Print out some solver info if desired:
if problem.options.verbose > 0
   disp('Transcription using GPOPS2');
end

% Copy the problem specification
setup = problem.options.gpops;

setup.bounds.phase.initialtime.lower = problem.bounds.initialTime.low';
setup.bounds.phase.initialtime.upper = problem.bounds.initialTime.upp';
setup.bounds.phase.finaltime.lower = problem.bounds.finalTime.low';
setup.bounds.phase.finaltime.upper = problem.bounds.finalTime.upp';

setup.bounds.phase.initialstate.lower = problem.bounds.initialState.low';
setup.bounds.phase.initialstate.upper = problem.bounds.initialState.upp';
setup.bounds.phase.finalstate.lower = problem.bounds.finalState.low';
setup.bounds.phase.finalstate.upper = problem.bounds.finalState.upp';
setup.bounds.phase.state.lower = problem.bounds.state.low';
setup.bounds.phase.state.upper = problem.bounds.state.upp';

setup.bounds.phase.control.lower = problem.bounds.control.low';
setup.bounds.phase.control.upper = problem.bounds.control.upp';

setup.guess.phase.time = problem.guess.time';
setup.guess.phase.state = problem.guess.state';
setup.guess.phase.control = problem.guess.control';

% Configure bounds on the path constraints
if ~isempty(problem.func.pathCst)
    if ~isfield(setup.bounds.phase, 'path')
        [cTest, ceqTest] =  problem.func.pathCst(...
            problem.guess.time,problem.guess.state,problem.guess.control);
        nc = size(cTest,1);
        nceq = size(ceqTest,1);
        setup.bounds.phase.path.lower = [-inf(1,nc), zeros(1,nceq)];
        setup.bounds.phase.path.upper = zeros(1,nc+nceq);
    end
end

% Configure bounds on the endpoint constraints
if ~isempty(problem.func.bndCst)
    if ~isfield(setup.bounds, 'eventgroup')
        t0 = problem.guess.time(1);  tF = problem.guess.time(end);
        x0 = problem.guess.state(:,1);  xF = problem.guess.state(:,end);
        [cTest, ceqTest] =  problem.func.bndCst(t0, x0,tF,xF);
        nc = size(cTest,1);
        nceq = size(ceqTest,1);
        setup.bounds.eventgroup.lower = [-inf(1,nc), zeros(1,nceq)];
        setup.bounds.eventgroup.upper = zeros(1,nc+nceq);
    end
end

F = problem.func;
setup.functions.continuous = @(input)( gpops_continuous(input,F.dynamics,F.pathObj,F.pathCst) );
setup.functions.endpoint = @(input)( gpops_endpoint(input,F.bndObj,F.bndCst) );

%%%% KEY LINE:  Solve the optimization problem with GPOPS II
output = gpops2(setup);

% Pack up the results:
soln.grid.time = output.result.solution.phase.time';
soln.grid.state = output.result.solution.phase.state';
soln.grid.control = output.result.solution.phase.control';

tSoln = output.result.interpsolution.phase.time';
xSoln = output.result.interpsolution.phase.state';
uSoln = output.result.interpsolution.phase.control';

soln.interp.state = @(t)( interp1(tSoln',xSoln',t','pchip',nan)' );
soln.interp.control = @(t)( interp1(tSoln',uSoln',t','pchip',nan)' );

soln.info.nlpTime = output.totaltime;
soln.info.objVal = output.result.objective;

soln.info.gpops.meshcounts = output.meshcounts;
soln.info.gpops.result.maxerror = output.result.maxerror;
soln.info.gpops.result.nlpinfo = output.result.nlpinfo;
soln.info.gpops.result.setup = output.result.setup;

soln.problem = problem;
soln.problem.options.nlpOpt = [];  % did not use the fmincon options

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function output = gpops_endpoint(input,bndObj,bndCst)
%
%  The endpoint function contains the boundary constraints and objective
%  functions for the trajectory optimization problem.
%

t0 = input.phase.initialtime;
tF = input.phase.finaltime;
x0 = input.phase.initialstate';
xF = input.phase.finalstate';

if isempty(bndObj)
    output.objective = input.phase.integral;
else
    output.objective = input.phase.integral + bndObj(t0,x0,tF,xF);
end

if ~isempty(bndCst)
    [c, ceq] = bndCst(t0,x0,tF,xF);
    output.eventgroup.event = [c;ceq]';
end

end



function output = gpops_continuous(input,dynamics,pathObj,pathCst)
%
% The continuous function contains the path objective, dynamics, and path
% constraint functions.
%

t = input.phase.time';
x = input.phase.state';
u = input.phase.control';
f = dynamics(t,x,u);
c = pathObj(t,x,u);

output.dynamics = f';
output.integrand = c';

if ~isempty(pathCst)
    [c, ceq] = pathCst(t,x,u);
    output.path = [c;ceq]';
end

end
