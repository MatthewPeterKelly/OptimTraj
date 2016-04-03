function problem = getDefaultOptions(problem)
% problem = getDefaultOptions(problem)
%
% This function fills in any blank entries in the problem.options struct.
% It is designed to be called from inside of trajOpt.m, and not by the
% user.
%

%%%% Top-level default options:
OPT.method = 'trapezoid';
OPT.verbose = 2;
OPT.defaultAccuracy = 'medium';


%%%% Basic setup

% ensure that options is not empty
if ~isfield(problem,'options')
    problem.options.method = OPT.method;
end
opt = problem.options;

% Loop over each options struct and fill in top-level options
for i=1:length(opt)
    if ~isfield(opt(i),'method')
        opt(i).method = OPT.method;
    elseif isempty(opt(i).method)
        opt(i).method = OPT.method;
    end
    if ~isfield(opt(i),'verbose')
        opt(i).verbose = OPT.verbose;
    elseif isempty(opt(i).verbose)
        opt(i).verbose = OPT.verbose;
    end
    if ~isfield(opt(i),'defaultAccuracy')
        opt(i).defaultAccuracy = OPT.defaultAccuracy;
    elseif isempty(opt(i).defaultAccuracy)
        opt(i).defaultAccuracy = OPT.defaultAccuracy;
    end
end

% Figure out basic problem size:
nState = size(problem.guess.state,1);
nControl = size(problem.guess.control,1);

% Loop over opt and fill in nlpOpt struct:
for i=1:length(opt)
    switch opt(i).verbose
        case 0
            NLP_display = 'notify';
        case 1
            NLP_display = 'final-detailed';
        case 2
            NLP_display = 'iter';
        case 3
            NLP_display = 'iter-detailed';
        otherwise
            error('Invalid value for options.verbose');
    end
    switch opt(i).defaultAccuracy
        case 'low'
            OPT.nlpOpt = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-4,...
                'MaxIter',100,...
                'MaxFunEvals',1000*(nState+nControl));
        case 'medium'
            OPT.nlpOpt = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-6,...
                'MaxIter',200,...
                'MaxFunEvals',2000*(nState+nControl));
        case 'high'
            OPT.nlpOpt = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-8,...
                'MaxIter',500,...
                'MaxFunEvals',5000*(nState+nControl));
        otherwise
            error('Invalid value for options.defaultAccuracy')
    end
    if isfield(opt(i),'nlpOpt')
        if isstruct(opt(i).nlpOpt) && ~isempty(opt(i).nlpOpt)
            names = fieldnames(opt(i).nlpOpt);
            for j=1:length(names)
                if ~isfield(OPT.nlpOpt,names{j})
                    disp(['WARNING: options.nlpOpt.' names{j} ' is not a valid option']);
                else
                    OPT.nlpOpt.(names{j}) = opt(i).nlpOpt.(names{j});
                end
            end
        end
    end
    opt(i).nlpOpt = OPT.nlpOpt;
end

% Check ChebFun dependency:
missingChebFun = false;
for i=1:length(opt)
if strcmp(opt(i).method,'chebyshev')
try
    chebpts(3);  %Test call to chebfun
catch ME %#ok<NASGU>
    missingChebFun = true;
    opt(i).method = 'trapezoid';  %Force default method
end
end
end
if missingChebFun
   warning('''chebyshev'' method requires the Chebfun toolbox');
   disp('   --> Install Chebfun toolbox:  (http://www.chebfun.org/)');
   disp('   --> Running with default method instead (''trapezoid'')');
end

% Fill in method-specific paramters:
for i=1:length(opt)
    OPT_method = opt(i).method;
    switch OPT_method
        case 'trapezoid'
            OPT.trapezoid = defaults_trapezoid(opt(i).defaultAccuracy);
        case 'hermiteSimpson'
            OPT.hermiteSimpson = defaults_hermiteSimpson(opt(i).defaultAccuracy);
        case 'chebyshev'
            OPT.chebyshev = defaults_chebyshev(opt(i).defaultAccuracy);
        case 'multiCheb'
            OPT.multiCheb = defaults_multiCheb(opt(i).defaultAccuracy);
        case 'rungeKutta'
            OPT.rungeKutta = defaults_rungeKutta(opt(i).defaultAccuracy);
        case 'gpops'
            OPT.gpops = defaults_gpops(opt(i).defaultAccuracy);
        otherwise
            error('Invalid value for options.method');
    end
    if isfield(opt(i),OPT_method)
        if isstruct(opt(i).(OPT_method)) && ~isempty(opt(i).(OPT_method))
            names = fieldnames(opt(i).(OPT_method));
            for j=1:length(names)
                if ~isfield(OPT.(OPT_method),names{j})
                    disp(['WARNING: options.' OPT_method '.' names{j} ' is not a valid option']);
                else
                    OPT.(OPT_method).(names{j}) = opt(i).(OPT_method).(names{j});
                end
            end
        end
    end
    opt(i).(OPT_method) = OPT.(OPT_method);
end

problem.options = opt;
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Method-specific parameters                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%



function OPT_trapezoid = defaults_trapezoid(accuracy)

switch accuracy
    case 'low'
        OPT_trapezoid.nGrid = 15;
    case 'medium'
        OPT_trapezoid.nGrid = 30;
    case 'high'
        OPT_trapezoid.nGrid = 60;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

end



function OPT_hermiteSimpson = defaults_hermiteSimpson(accuracy)

switch accuracy
    case 'low'
        OPT_hermiteSimpson.nSegment = 10;
    case 'medium'
        OPT_hermiteSimpson.nSegment = 20;
    case 'high'
        OPT_hermiteSimpson.nSegment = 40;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

end



function OPT_chebyshev = defaults_chebyshev(accuracy)

switch accuracy
    case 'low'
        OPT_chebyshev.nColPts = 9;
    case 'medium'
        OPT_chebyshev.nColPts = 13;
    case 'high'
        OPT_chebyshev.nColPts = 23;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

end


function OPT_multiCheb = defaults_multiCheb(accuracy)

switch accuracy
    case 'low'
        OPT_multiCheb.nColPts = 6;
        OPT_multiCheb.nSegment = 3;
    case 'medium'
        OPT_multiCheb.nColPts = 8;
        OPT_multiCheb.nSegment = 6;
    case 'high'
        OPT_multiCheb.nColPts = 8;
        OPT_multiCheb.nSegment = 12;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

end


function OPT_rungeKutta = defaults_rungeKutta(accuracy)

switch accuracy
    case 'low'
        OPT_rungeKutta.nSegment = 10;
        OPT_rungeKutta.nSubStep = 2;
    case 'medium'
        OPT_rungeKutta.nSegment = 20;
        OPT_rungeKutta.nSubStep = 2;
    case 'high'
        OPT_rungeKutta.nSegment = 20;
        OPT_rungeKutta.nSubStep = 4;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

OPT_rungeKutta.AdaptiveDerivativeCheck = 'off';

end




function OPT_gpops = defaults_gpops(accuracy)

OPT_gpops.bounds.phase.integral.lower = -inf;
OPT_gpops.bounds.phase.integral.upper = inf;
OPT_gpops.guess.phase.integral = 0;

OPT_gpops.name = 'TrajOpt_GPOPS';
OPT_gpops.auxdata = [];
OPT_gpops.nlp.solver = 'ipopt'; % {'ipopt','snopt'}
OPT_gpops.derivatives.dependencies = 'full';  %�full�, �sparse� or �sparseNaN�
OPT_gpops.derivatives.supplier = 'sparseCD'; %'sparseCD';  %'adigator'
OPT_gpops.derivatives.derivativelevel = 'first'; %'second';
OPT_gpops.mesh.method = 'hp-PattersonRao';
OPT_gpops.method = 'RPM-Integration';
OPT_gpops.mesh.phase.colpoints = 10*ones(1,10);
OPT_gpops.mesh.phase.fraction = ones(1,10)/10;
OPT_gpops.scales.method = 'none'; % { 'none' , automatic-hybridUpdate' , 'automatic-bounds';

switch accuracy
    case 'low'
        OPT_gpops.mesh.tolerance = 1e-2;
        OPT_gpops.mesh.maxiterations = 0;
    case 'medium'
        OPT_gpops.mesh.tolerance = 1e-3;
        OPT_gpops.mesh.maxiterations = 1;
        
    case 'high'
        OPT_gpops.mesh.tolerance = 1e-4;
        OPT_gpops.mesh.maxiterations = 3;
    otherwise
        error('Invalid value for options.defaultAccuracy')
end

end
