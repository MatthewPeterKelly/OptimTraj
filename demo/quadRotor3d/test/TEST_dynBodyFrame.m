disp('Running: TEST_dynBodyFrame.m')
% define parameter struct
p.g = -9.81 ; 
p.rho = 1.1225 ; 
p.m = 5 ; 
p.I = [0.625 0 0; 0 0.625 0; 0 0 1.25] ; % inertia tensor
p.cg = [0 0 0] ; % (m) location of center of gravity
qRP.d_prop = 0.305*ones(4,1) ; % propeller diameter (m)
qRP.maxThrust = 25*ones(4,1) ; % thrust at 100% throttle (N)
qRP.maxRPM = 10000*ones(4,1) ; % RPM at 100% throttle (RPM)
qRP.maxTorque = ones(4,1) ;  % torque at 100% throttle (Nm)
qRP.thrustLocations = [0.5 0 0; 0 0.5 0; -0.5 0 0; 0 -0.5 0]; % motor locations (each row one motor in coords: [port, nose, top] 
qRP.thrustAxes = repmat([0 0 1],4,1) ; % thrust axes of each motor in coords port, nose, top.
qRP.isSpinDirectionCCW = [1; 0; 1; 0] ; % bool to reverse motor spin direction around 'thrustAxes'.

% Call function that creates the propulsion plant model
[p.propulsion] = definePropulsionModel(qRP); clear variable qRP; 

%% Test 0 - No input force => no acceleration
disp('Test 0 - no input force')
u = zeros(4,1) ; 
omega = zeros(3,1) ; 
[ds] = dynBodyFrame(omega, u, p) ;

%% Test 1 - No input force => no acceleration
disp('Test 1 - zero throttle ''wide vector'' ') 
u = zeros(4,10) ; 
omega = zeros(3,10) ; 
[ds] = dynBodyFrame(omega,u, p) ;

%% Test 2 - All motors at max throttle => acceleration in positive Z direction [i.e. ds(3)]
disp('Test 2 - max throttle') 
u = ones(4,10) ; 
omega = zeros(3,10) ; 
[ds] = dynBodyFrame(omega,u, p) ; 

%% Test 3 - "wide" input vector.
disp('Test 3 - ramp throttle ''wide'' ') 
u = repmat(linspace(0,1,100),4,1); 
omega = rand(3,100) ; 
[ds] = dynBodyFrame(omega,u, p) ;

%% Test 4 - only motor 1 (port side motor {at the moment}) 
% should result in (pitch accel = 0, roll accel < 0, yaw accel < 0)
disp('Test 4 - port motor only ; ramp') 
u = zeros(4,100) ;
u(1,:) = linspace(0,1,100) ; 
omega = zeros(3,100) ; 
[ds] = dynBodyFrame(omega,u, p) ;

%% Test 5 - only motor 2 (nose motor {at the moment}) 
% should result in (pitch accel > 0, roll accel = 0, yaw accel > 0)
disp('Test 5 - nose motor only ; ramp') 
u = zeros(4,100) ;
u(2,:) = linspace(0,1,100) ; 
omega = zeros(3,100) ; 
[ds] = dynBodyFrame(omega,u, p) ;

%% 
disp('TEST_dynBodyFrame concluded without errors')