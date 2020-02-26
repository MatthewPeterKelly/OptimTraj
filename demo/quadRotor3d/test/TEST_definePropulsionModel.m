% defines an aircraft propulsion model
disp('Running: TEST_definePropulsionModel.m') 

% Define example parameters
qRP.d_prop = 0.305*ones(4,1) ; % propeller diameter (m)
qRP.maxThrust = 25*ones(4,1) ; % thrust at 100% throttle (N)
qRP.maxRPM = 10000*ones(4,1) ; % RPM at 100% throttle (RPM)
qRP.maxTorque = ones(4,1) ;  % torque at 100% throttle (Nm)
qRP.thrustLocations = [0.5 0 0; 0 0.5 0; -0.5 0 0; 0 -0.5 0]; % motor locations (each row one motor in coords: [port, nose, top] 
qRP.thrustAxes = repmat([0 0 1],4,1) ; % thrust axes of each motor in coords port, nose, top.
qRP.isSpinDirectionCCW = [1; 0; 1; 0] ; % bool to reverse motor spin direction around 'thrustAxes'.

% Call function that creates the propulsion plant model (without plotting)
[quadrotorPropulsionModel_1] = definePropulsionModel(qRP) ;

% Call function that creates the propulsion plant model (with plotting)
plotflag = 1; 
[quadrotorPropulsionModel_2] = definePropulsionModel(qRP,plotflag) ;

%%
disp('TEST_definePropulsionModel.m ran without error') 