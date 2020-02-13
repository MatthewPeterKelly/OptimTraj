% defines an aircraft propulsion model

% Define example parameters
qRP.d_prop = 0.305 ; % propeller diameter (m)
qRP.maxThrust = 25 ; % thrust at 100% throttle (N)
qRP.maxRPM = 10000 ; % RPM at 100% throttle (RPM)
qRP.maxTorque = 1 ;  % torque at 100% throttle (Nm)
qRP.thrustLocations = [0.5 0 0; 0 0.5 0; -0.5 0 0; 0 -0.5 0]; % motor locations (each row one motor in coords: [port, nose, top] 
qRP.thrustAxes = repmat([0 0 1],4,1) ; % thrust axes of each motor in coords port, nose, top.
qRP.isSpinDirectionCCW = [1; 0; 1; 0] ; % bool to reverse motor spin direction around 'thrustAxes'.

% Call function that creates the propulsion plant model
[quadrotorPropulsionModel] = definePropulsionModel(qRP)

% Plot newly created model
showPropulsionModel(quadrotorPropulsionModel)