disp('Running: TEST_dynQuadRotor3d.m')
s = zeros(12,1) ; % [x, y, z, pitch, roll, yaw, dx, dy, dz, dpitch, droll, dyaw] 
u = zeros(4,1) ; 

% Coordinate Frame Note:
% World frame is XYZ = [East, North, Up]. 
% Body frame is XYZ = [port, nose, top], so attitude Euler angles are in
% the order: [pitch, roll, yaw]


% Define environmental and plant model params
% Enviromental params
p.g = -9.81 ; % World Coords is XYZ = [East, North, Up], i.e. gravity is a negative number
p.rho = 1.225 ; 

% Inertial params
p.m = 5 ; 
p.I = [0.625 0 0; 0 0.625 0; 0 0 1.25] ; % inertia tensor
p.cg = [0 0 0] ; % (m) location of center of gravity

% Propulsion params
define_propulsion_model
p.propulsion = propulsion ; clear variable propulsion ; 

%% Single Timestep
% Throttle = 0 %
disp('Test 0 - 0% throttle')
[ds] = dynQuadRotor3d(s, u, p)

%% Single timestep - full throttle
disp('Test 1 - 100% throttle')
s = zeros(12,1) ;
u = ones(4,1) ; 
[ds] = dynQuadRotor3d(s, u, p)

%% multiple timesteps
disp('Test 2 - 0% throttle')
u = zeros(4,10) ; 
s = zeros(12,10) ; 
[ds] = dynQuadRotor3d(s, u, p)

%% multiple timesteps
disp('Test 3 - 0% throttle')
u = ones(4,10) ; 
s = zeros(12,10) ;
s(4,:) = linspace(0,1,10) ; 
[ds] = dynQuadRotor3d(s, u, p)

%% Effect of aircraft orientation
% positive pitch; should result in an acceleration 'south' (i.e. negative
% north => ds(8) < 0 ) ; 
disp('Test 4 - 100% throttle, pitch = 5 deg')
deflect = 5 ; 
s = zeros(12,10) ;
s(4,:) = deg2rad(deflect) ;
u = ones(4,10) ;
[ds] = dynQuadRotor3d(s, u, p)

%% negative pitch; should result in an acceleration 'north'
disp('Test 3 - 100% throttle, pitch = -5 deg')
deflect = -5 ; 
s = zeros(12,10) ;
s(4,:) = deg2rad(deflect) 
u = ones(4,10)  
[ds] = dynQuadRotor3d(s, u, p)

%% positive roll; should result in an acceleration 'east', i.e. ds(7) > 0
disp('Test 4 - 100% throttle, roll = 5 deg')
deflect = 5 ; 
s = zeros(12,10) ;
s(5,:) = deg2rad(deflect) 
u = ones(4,10)  
[ds] = dynQuadRotor3d(s, u, p)

%% negative roll; should result in an acceleration 'west', i.e. negative element 7
disp('Test 5 - 100% throttle, roll = 5 deg')
deflect = -5 ; 
s = zeros(12,1) ;
s(5) = deg2rad(deflect) 
u = ones(4,1)  
[ds] = dynQuadRotor3d(s, u, p)

%% positive yaw; should not affect acceleration direction
disp('Test 6 - 100% throttle, roll = 5 deg')
deflect = 5 ; 
s = zeros(12,1) ;
s(6) = deg2rad(deflect) ;
u = ones(4,1)   ;
[ds] = dynQuadRotor3d(s, u, p) ;

%% negative yaw; should not affect acceleration direction
disp('Test 7 - 100% throttle, roll = 5 deg')
deflect = -5 ; 
s = zeros(12,1) ;
s(6) = deg2rad(deflect) ; 
u = ones(4,1)  ;
[ds] = dynQuadRotor3d(s, u, p) 

%% wide vector
s = zeros(12,100) ;
u = repmat(linspace(0,1,100),4,1) ; 
[ds] = dynQuadRotor3d(s, u, p)
