disp('Running: TEST_computePropOpPoint.m')

%% Single operating points
% define some basic parameters
rho = 1.225 ; 
d_prop = 0.305 ; 
C_t = 0.0849 ; 
C_q = 0.0111 ; 

%% RPM = 0 
disp('Test 1 - RPM = 0')
RPM = 0 ; 
[thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q) ;

%% RPM = 1000 
disp('Test 2 - RPM = 1000')
RPM = 1000 ; 
[thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q) ;

%% Multiple operating points
disp('Test 3 - RPM = [1000 2000 3000]')
RPM = [1000 2000 3000] ; 
[thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q) ;

%%
disp('TEST_computePropOpPoint.m ran without error') 