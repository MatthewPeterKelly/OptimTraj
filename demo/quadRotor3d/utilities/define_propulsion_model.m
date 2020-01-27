% Defines a propulsion parameter plant model.

% Assumes same prop-motor combination at all locations.

propulsion = struct() ; % initialize struct
id = 0 ; % use counter to make code reuse easier

%% Propulsion system test data (at 100% throttle)
% Change this for a different propulsion system
d_prop = 0.305 ; % propeller diameter (m)
maxThrust = 25 ; % thrust at 100% throttle (N)
maxRPM = 10000 ; % RPM at 100% throttle (RPM)
maxTorque = 1 ;  % torque at 100% throttle (Nm)
rho = 1.225 ;    % air density during data collection (kg/m^3) 

% Compute propulsion system coefficients
% See: https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
C_t = maxThrust / (rho * (maxRPM/60)^2 * d_prop^4) ; 
C_q = maxTorque / (rho * (maxRPM/60)^2 * d_prop^5) ;

%% Assign to struct
% propeller 1
id = id + 1 ; 
propulsion(id).thrustAxis = [0 0 1] ;      % [port, nose, top] % nose faces north when body and world axis are aligned (world is "East North Up")
propulsion(id).thrustLocation = [0.5 0 0] ; 
propulsion(id).isSpinDirectionCCW = 1 ; 
propulsion(id).C_t = C_t ; % initialize
propulsion(id).C_q = C_q ; % initialize
propulsion(id).maxRPM = maxRPM ; 
propulsion(id).maxTorque = maxTorque ; 
propulsion(id).d_prop = d_prop ; 

% propeller 2
id = id + 1 ; 
propulsion(id).thrustAxis = [0 0 1] ; 
propulsion(id).thrustLocation = [0 0.5 0] ; 
propulsion(id).isSpinDirectionCCW = 0 ; 
propulsion(id).C_t = C_t ; % initialize
propulsion(id).C_q = C_q ; % initialize
propulsion(id).maxRPM = maxRPM ; 
propulsion(id).maxTorque = maxTorque ; 
propulsion(id).d_prop = d_prop ; 

% propeller 3
id = id + 1 ; 
propulsion(id).thrustAxis = [0 0 1] ; 
propulsion(id).thrustLocation = [-0.5 0 0] ; 
propulsion(id).isSpinDirectionCCW = 1 ; 
propulsion(id).C_t = C_t ; % initialize
propulsion(id).C_q = C_q ; % initialize
propulsion(id).maxRPM = maxRPM ; 
propulsion(id).maxTorque = maxTorque ; 
propulsion(id).d_prop = d_prop ; 

% propeller 4
id = id + 1 ; 
propulsion(id).thrustAxis = [0 0 1] ; 
propulsion(id).thrustLocation = [0 -0.5 0] ; 
propulsion(id).isSpinDirectionCCW = 0 ; 
propulsion(id).C_t = C_t ; % initialize
propulsion(id).C_q = C_q ; % initialize
propulsion(id).maxRPM = maxRPM ; 
propulsion(id).maxTorque = maxTorque ; 
propulsion(id).d_prop = d_prop ; 

%% cleanup
clear variable C_q C_t d_prop id maxRPM maxThrust maxTorque