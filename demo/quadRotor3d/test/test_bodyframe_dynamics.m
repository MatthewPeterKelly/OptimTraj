disp('Running: test_bodyframe_dynamics.m')
% define parameter struct
p.g = -9.81 ; 
p.rho = 1.1225 ; 
p.m = 5 ; 
p.I = [0.625 0 0; 0 0.625 0; 0 0 1.25] ; % inertia tensor
p.cg = [0 0 0] ; % (m) location of center of gravity
define_propulsion_model
p.propulsion = propulsion ; clear variable propulsion ; 

% Test 0 - No input force => no acceleration
disp('Test 0') 
u = zeros(4,1) ; 
[ds] = bodyframe_dynamics(u, p)

% Test 1 - All motors at max throttle => acceleration in positive Z direction [i.e. ds(3)]
disp('Test 1') 
u = ones(4,1) ; 
[ds] = bodyframe_dynamics(u, p)

%% Test 2 - "wide" input vector.
disp('Test 2') 
u = repmat(linspace(0,1,100),4,1); 
[ds] = bodyframe_dynamics(u, p)
figure
plot(ds(3,:),'o') ; grid on ; hold on; 

%% Test 3 - only motor 1 (port side motor {at the moment}) 
% should result in (pitch accel = 0, roll accel < 0, yaw accel < 0)
disp('Test 3') 
u = zeros(4,100) ;
u(1,:) = linspace(0,1,100) ; 
[ds] = bodyframe_dynamics(u, p)

%% Test 4 - only motor 2 (nose motor {at the moment}) 
% should result in (pitch accel > 0, roll accel = 0, yaw accel > 0)
disp('Test 4') 
u = zeros(4,100) ;
u(2,:) = linspace(0,1,100) ; 
[ds] = bodyframe_dynamics(u, p)