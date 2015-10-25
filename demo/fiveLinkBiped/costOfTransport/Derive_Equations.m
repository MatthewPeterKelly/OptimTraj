function Derive_Equations()
%%%% Derive Equations - Five Link Biped Model %%%%
%
% This function derives the equations of motion, as well as some other useful
% equations (kinematics, contact forces, ...) for the five-link biped
% model.
%
% This version of the code includes a few more complicated features for
% dealing with difficult cost functions. In particular, it adds 10 slack
% variables to compute the abs(power) term in the cost function, and the
% primary control is the derivative of torque, rather than torque itself.
% This allows for regularization by the derivative of the input.
%
%
% Nomenclature:
%
% - There are five links, which will be numbered starting with "1" for the
% stance leg tibia, increasing as the links are father from the base joint,
% and ending with "5" for the swing leg tibia.
%   1 - stance leg tibia (lower leg)
%   2 - stance leg femur  (upper leg)
%   3 - torso
%   4 - swing leg femur
%   5 - swing leg tibia
%
% - This script uses absolute angles, which are represented with "q". All
% angles use positive convention, with the zero angle corresponding to a
% vertically aligned link configuration. [q] = [0] has the torso balanced
% upright, with both legs fully extended straight below it.
%
% - Derivatives with respect to time are notated by prepending a "d". For
% example the rate of change in an absolute angle is "dq" and angular
% acceleration would be "ddq"
%
% - Joint positions are given with "P", center of mass positions are "G"
%

clc; clear;
disp('Creating variables and derivatives...')

%%%% Absolute orientation (angle) of each link
q1 = sym('q1', 'real');
q2 = sym('q2','real');
q3 = sym('q3','real');
q4 = sym('q4','real');
q5 = sym('q5','real');

%%%% Absolute angular rate of each link
dq1 = sym('dq1','real');
dq2 = sym('dq2','real');
dq3 = sym('dq3','real');
dq4 = sym('dq4','real');
dq5 = sym('dq5','real');

%%%% Absolute angular acceleration of each linke
ddq1 = sym('ddq1','real');
ddq2 = sym('ddq2','real');
ddq3 = sym('ddq3','real');
ddq4 = sym('ddq4','real');
ddq5 = sym('ddq5','real');

%%%% Torques at each joint
u1 = sym('u1','real');  %Stance foot
u2 = sym('u2','real');   %Stance knee
u3 = sym('u3','real');   %Stance hip
u4 = sym('u4','real');   %Swing hip
u5 = sym('u5','real');   %Swing knee

%%%% Torques rate at each joint
du1 = sym('du1','real');  %Stance foot
du2 = sym('du2','real');   %Stance knee
du3 = sym('du3','real');   %Stance hip
du4 = sym('du4','real');   %Swing hip
du5 = sym('du5','real');   %Swing knee

%%%% Slack variables -- negative component of power
sn1 = sym('sn1','real');  %Stance foot
sn2 = sym('sn2','real');   %Stance knee
sn3 = sym('sn3','real');   %Stance hip
sn4 = sym('sn4','real');   %Swing hip
sn5 = sym('sn5','real');   %Swing knee

%%%% Slack variables -- positive component of power
sp1 = sym('sp1','real');  %Stance foot
sp2 = sym('sp2','real');   %Stance knee
sp3 = sym('sp3','real');   %Stance hip
sp4 = sym('sp4','real');   %Swing hip
sp5 = sym('sp5','real');   %Swing knee

%%%% Mass of each link
m1 = sym('m1','real');
m2 = sym('m2','real');
m3 = sym('m3','real');
m4 = sym('m4','real');
m5 = sym('m5','real');

%%%% Distance between parent joint and link center of mass
c1 = sym('c1','real');
c2 = sym('c2','real');
c3 = sym('c3','real');
c4 = sym('c4','real');
c5 = sym('c5','real');

%%%% Length of each link
l1 = sym('l1','real');
l2 = sym('l2','real');
l3 = sym('l3','real');
l4 = sym('l4','real');
l5 = sym('l5','real');

%%%% Moment of inertia of each link about its own center of mass
I1 = sym('I1','real');
I2 = sym('I2','real');
I3 = sym('I3','real');
I4 = sym('I4','real');
I5 = sym('I5','real');

g = sym('g','real'); % Gravity
Fx = sym('Fx','real');   %Horizontal contact force at stance foot
Fy = sym('Fy','real');   %Vertical contact force at stance foot
empty = sym('empty','real');   %Used for vectorization, user should pass a vector of zeros
t = sym('t','real');  %dummy continuous time

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Set up coordinate system and unit vectors                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

i = sym([1;0]);   %Horizontal axis
j = sym([0;1]);   %Vertical axis

e1 = cos(q1)*(j) + sin(q1)*(-i);  %unit vector from P0 -> P1, (contact point to stance knee)
e2 = cos(q2)*(j) + sin(q2)*(-i);  %unit vector from P1 -> P2, (stance knee to hip)
e3 = cos(q3)*(j) + sin(q3)*(-i);  %unit vector from P2 -> P3, (hip to shoulders);
e4 = -cos(q4)*(j) - sin(q4)*(-i);  %unit vector from P2 -> P4, (hip to swing knee);
e5 = -cos(q5)*(j) - sin(q5)*(-i);  %unit vector from P4 -> P5, (swing knee to swing foot);

P0 = 0*i + 0*j;   %stance foot = Contact point = origin
P1 = P0 + l1*e1;  %stance knee
P2 = P1 + l2*e2;  %hip
P3 = P2 + l3*e3;  %shoulders
P4 = P2 + l4*e4;  %swing knee
P5 = P4 + l5*e5;  %swing foot

G1 = P1 - c1*e1;  % CoM stance leg tibia
G2 = P2 - c2*e2;  % CoM stance leg febur
G3 = P3 - c3*e3;  % CoM torso
G4 = P2 + c4*e4;  % CoM swing leg femur
G5 = P4 + c5*e5;  % CoM swing leg tibia
G = (m1*G1 + m2*G2 + m3*G3 + m4*G4 + m5*G5)/(m1+m2+m3+m4+m5);  %Center of mass for entire robot

%%%% Define a function for doing '2d' cross product: dot(a x b, k)
cross2d = @(a,b)(a(1)*b(2) - a(2)*b(1));

%%%% Weight of each link:
w1 = -m1*g*j;
w2 = -m2*g*j;
w3 = -m3*g*j;
w4 = -m4*g*j;
w5 = -m5*g*j;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                             Derivatives                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

q = [q1;q2;q3;q4;q5];
dq = [dq1;dq2;dq3;dq4;dq5];
ddq = [ddq1;ddq2;ddq3;ddq4;ddq5];
u = [u1;u2;u3;u4;u5];
du = [du1;du2;du3;du4;du5];
sn = [sn1;sn2;sn3;sn4;sn5];
sp = [sp1;sp2;sp3;sp4;sp5];
z = [t;q;dq;u;du;sn;sp];   % time-varying vector of inputs

% Neat trick to compute derivatives using the chain rule
derivative = @(in)( jacobian(in,[q;dq;u])*[dq;ddq;du] );

% Velocity of the swing foot (used for step constraints)
dP5 = derivative(P5);

% Compute derivatives for the CoM of each link:
dG1 = derivative(G1);  ddG1 = derivative(dG1);
dG2 = derivative(G2);  ddG2 = derivative(dG2);
dG3 = derivative(G3);  ddG3 = derivative(dG3);
dG4 = derivative(G4);  ddG4 = derivative(dG4);
dG5 = derivative(G5);  ddG5 = derivative(dG5);
dG = derivative(G);  ddG = derivative(dG);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                             Calculations:                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

singleStanceDynamics();
objectiveFunctions();
heelStrikeDynamics();

mechanicalEnergy();
contactForces();
kinematics();

disp('Done!');

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Single-Stance Dynamics                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% I solve the dynamics here by carefully selecting angular momentum balance
% equations about each joint, working my way out the kinematic tree from
% the root.

    function singleStanceDynamics()
        disp('Deriving single stance dynamics...')
        
        %%%% AMB - entire system @ P0
        eqnTorque0 = ...
            cross2d(G1-P0,w1) + ...
            cross2d(G2-P0,w2) + ...
            cross2d(G3-P0,w3) + ...
            cross2d(G4-P0,w4) + ...
            cross2d(G5-P0,w5) + ...
            u1;
        
        eqnInertia0 = ...
            cross2d(G1-P0,m1*ddG1) + ddq1*I1 + ...
            cross2d(G2-P0,m2*ddG2) + ddq2*I2 + ...
            cross2d(G3-P0,m3*ddG3) + ddq3*I3 + ...
            cross2d(G4-P0,m4*ddG4) + ddq4*I4 + ...
            cross2d(G5-P0,m5*ddG5) + ddq5*I5;
        
        %%%% AMB - swing leg, torso, stance femer  @ stance knee
        eqnTorque1 = ...
            cross2d(G2-P1,w2) + ...
            cross2d(G3-P1,w3) + ...
            cross2d(G4-P1,w4) + ...
            cross2d(G5-P1,w5) + ...
            u2;
        
        eqnInertia1 = ...
            cross2d(G2-P1,m2*ddG2) + ddq2*I2  + ...
            cross2d(G3-P1,m3*ddG3) + ddq3*I3  + ...
            cross2d(G4-P1,m4*ddG4) + ddq4*I4  + ...
            cross2d(G5-P1,m5*ddG5) + ddq5*I5 ;
        
        %%%% AMB - swing leg, torso @ hip
        eqnTorque2 = ...
            cross2d(G3-P2,w3) + ...
            cross2d(G4-P2,w4) + ...
            cross2d(G5-P2,w5) + ...
            u3;
        
        eqnInertia2 = ...
            cross2d(G3-P2,m3*ddG3) + ddq3*I3  + ...
            cross2d(G4-P2,m4*ddG4) + ddq4*I4  + ...
            cross2d(G5-P2,m5*ddG5) + ddq5*I5 ;
        
        %%%% AMB - swing leg @ hip
        eqnTorque3 = ...
            cross2d(G4-P2,w4) + ...
            cross2d(G5-P2,w5) + ...
            u4;
        
        eqnInertia3 = ...
            cross2d(G4-P2,m4*ddG4) + ddq4*I4  + ...
            cross2d(G5-P2,m5*ddG5) + ddq5*I5 ;
        
        %%%% AMB - swing tibia % swing knee
        eqnTorque4 = ...
            cross2d(G5-P4,w5) + ...
            u5;
        
        eqnInertia4 = ...
            cross2d(G5-P4,m5*ddG5) + ddq5*I5 ;
        
        %%%% Collect and solve equations:
        eqns = [...
            eqnTorque0 - eqnInertia0;
            eqnTorque1 - eqnInertia1;
            eqnTorque2 - eqnInertia2;
            eqnTorque3 - eqnInertia3;
            eqnTorque4 - eqnInertia4];
        
        [MM, FF] = equationsToMatrix(eqns,ddq);  % ddq = MM\ff;
        
        %%%% Compute gradients:
        [m, mi, mz, mzi, mzd] = computeGradients(MM,z,empty);
        [f, fi, fz, fzi, fzd] = computeGradients(FF,z,empty);
        
        % Write function file:
        matlabFunction(m, mi, f, fi,...   %dynamics
            mz, mzi, mzd, fz, fzi, fzd,...  %gradients
            'file','autoGen_dynSs.m',...
            'vars',{...
            'q1','q2','q3','q4','q5',...
            'dq1','dq2','dq3','dq4','dq5',...
            'u1','u2','u3','u4','u5',...
            'm1','m2','m3','m4','m5',...
            'I1','I2','I3','I4','I5',...
            'l1','l2','l3','l4',...
            'c1','c2','c3','c4','c5',...
            'g','empty'});
        
    end




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Objective Functions                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    function objectiveFunctions()
        
        % Joint rates:
        v1 = dq1;   % joint rate 1
        v2 = dq2-dq1;   % joint rate 2
        v3 = dq3-dq2; % joint rate 3
        v4 = dq4-dq3;  % joint rate 4
        v5 = dq5-dq4;  % joint rate 5
        
        % Compute the power used by each joint
        pow1 = v1*u1;  %Power used by joint 1
        pow2 = v2*u2;  %Power used by joint 2
        pow3 = v3*u3;  %Power used by joint 3
        pow4 = v4*u4;  %Power used by joint 4
        pow5 = v5*u5;  %Power used by joint 5
        
        % Constraint on the slack variables:
        slackCst = [...
            pow1 - (sp1 - sn1);
            pow2 - (sp2 - sn2);
            pow3 - (sp3 - sn3);
            pow4 - (sp4 - sn4);
            pow5 - (sp5 - sn5)];
        
        % Gradients of the constraint on slack variables:
        [c, ~, cz, czi, ~] = computeGradients(slackCst,z,empty);
        
        matlabFunction(c,cz,czi,...
            'file','autoGen_cst_costOfTransport.m',...
            'vars',{...
            'dq1','dq2','dq3','dq4','dq5',...
            'u1','u2','u3','u4','u5',...
            'sn1','sn2','sn3','sn4','sn5',...
            'sp1','sp2','sp3','sp4','sp5','empty'});
        
        % abs(power) using slack variables:
        gammaNeg = sym('gammaNeg','real');  %Torque-squared smoothing parameter
        gammaPos = sym('gammaPos','real');  %Torque-squared smoothing parameter
        absPower = gammaNeg*(sn1 + sn2 + sn3 + sn4 + sn5) + ...
            gammaPos*(sp1 + sp2 + sp3 + sp4 + sp5);
        
        % Cost of Transport:
        weight = (m1+m2+m3+m4+m5)*g;
        stepLength = sym('stepLength','real');
        alpha = sym('alpha','real');  %Torque-squared smoothing parameter
        beta = sym('beta','real');  %Torque-rate squared smoothing
        F = absPower/(weight*stepLength) + ...
            alpha*(u1^2 + u2^2 + u3^2 + u4^2 + u5^2) + ...
            beta*(du1^2 + du2^2 + du3^2 + du4^2 + du5^2);
        [f, ~, fz, fzi, ~]  = computeGradients(F,z,empty);
        
        matlabFunction(f,fz,fzi,...
            'file','autoGen_obj_costOfTransport.m',...
            'vars',{...
            'm1','m2','m3','m4','m5',...
            'u1','u2','u3','u4','u5',...
            'du1','du2','du3','du4','du5',...
            'sn1','sn2','sn3','sn4','sn5', ...
            'sp1','sp2','sp3','sp4','sp5',...
            'g','stepLength','gammaNeg','gammaPos','alpha','beta','empty'});
        
        % Swing foot height:
        stepHeight = sym('stepHeight','real');
        yFoot = P5(2);
        xFoot = P5(1);
        yMin = stepHeight*(1 - (xFoot/stepLength)^2);
        yCst = yMin - yFoot;  %Must be negative
        [y, ~, yz, yzi, ~]  = computeGradients(yCst,z,empty);
        matlabFunction(y,yz,yzi,...
                    'file','autoGen_cst_swingFootHeight.m',...
            'vars',{...
            'q1','q2','q4','q5',...
            'l1','l2','l4','l5'...
            'stepLength','stepHeight'});
    end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Heel-Strike Dynamics                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    function heelStrikeDynamics()
        disp('Deriving heel-strike dynamics...')
        
        %%%% Notes:
        % xF - heelStrike(xI) --> constraint --> 0
        % xF - collision(footSwap(xI));
        %
        
        % Angles before heel-strike:
        q1m = sym('q1m','real');
        q2m = sym('q2m','real');
        q3m = sym('q3m','real');
        q4m = sym('q4m','real');
        q5m = sym('q5m','real');
        qm = [q1m;q2m;q3m;q4m;q5m];
        
        % Angles after heel-strike
        q1p = sym('q1p','real');
        q2p = sym('q2p','real');
        q3p = sym('q3p','real');
        q4p = sym('q4p','real');
        q5p = sym('q5p','real');
        qp = [q1p;q2p;q3p;q4p;q5p];
        
        % Angular rates before heel-strike:
        dq1m = sym('dq1m','real');
        dq2m = sym('dq2m','real');
        dq3m = sym('dq3m','real');
        dq4m = sym('dq4m','real');
        dq5m = sym('dq5m','real');
        dqm = [dq1m;dq2m;dq3m;dq4m;dq5m];
        
        % Angular rates after heel-strike
        dq1p = sym('dq1p','real');
        dq2p = sym('dq2p','real');
        dq3p = sym('dq3p','real');
        dq4p = sym('dq4p','real');
        dq5p = sym('dq5p','real');
        dqp = [dq1p;dq2p;dq3p;dq4p;dq5p];
        
        % torque before heel-strike:
        u1m = sym('u1m','real');
        u2m = sym('u2m','real');
        u3m = sym('u3m','real');
        u4m = sym('u4m','real');
        u5m = sym('u5m','real');
        um = [u1m;u2m;u3m;u4m;u5m];
        
        % torque after heel-strike
        u1p = sym('u1p','real');
        u2p = sym('u2p','real');
        u3p = sym('u3p','real');
        u4p = sym('u4p','real');
        u5p = sym('u5p','real');
        up = [u1p;u2p;u3p;u4p;u5p];
        
        % Compute kinematics before heel-strike:
        inVars = {'q1','q2','q3','q4','q5','dq1','dq2','dq3','dq4','dq5'};
        outVarsM = {'q1m','q2m','q3m','q4m','q5m','dq1m','dq2m','dq3m','dq4m','dq5m'};
        %         P0m = subs(P0,inVars,outVarsM);
        P1m = subs(P1,inVars,outVarsM);
        P2m = subs(P2,inVars,outVarsM);
        %         P3m = subs(P3,inVars,outVarsM);
        P4m = subs(P4,inVars,outVarsM);
        P5m = subs(P5,inVars,outVarsM);
        dP5m = subs(dP5,inVars,outVarsM);
        G1m = subs(G1,inVars,outVarsM);
        G2m = subs(G2,inVars,outVarsM);
        G3m = subs(G3,inVars,outVarsM);
        G4m = subs(G4,inVars,outVarsM);
        G5m = subs(G5,inVars,outVarsM);
        dG1m = subs(dG1,inVars,outVarsM);
        dG2m = subs(dG2,inVars,outVarsM);
        dG3m = subs(dG3,inVars,outVarsM);
        dG4m = subs(dG4,inVars,outVarsM);
        dG5m = subs(dG5,inVars,outVarsM);
        
        % Compute kinematics after heel-strike:
        outVarsP = {'q1p','q2p','q3p','q4p','q5p','dq1p','dq2p','dq3p','dq4p','dq5p'};
        P0p = subs(P0,inVars,outVarsP);
        P1p = subs(P1,inVars,outVarsP);
        P2p = subs(P2,inVars,outVarsP);
        %         P3p = subs(P3,inVars,outVarsP);
        P4p = subs(P4,inVars,outVarsP);
        %         P5p = subs(P5,inVars,outVarsP);
        dP5p = subs(dP5,inVars,outVarsP);
        G1p = subs(G1,inVars,outVarsP);
        G2p = subs(G2,inVars,outVarsP);
        G3p = subs(G3,inVars,outVarsP);
        G4p = subs(G4,inVars,outVarsP);
        G5p = subs(G5,inVars,outVarsP);
        dG1p = subs(dG1,inVars,outVarsP);
        dG2p = subs(dG2,inVars,outVarsP);
        dG3p = subs(dG3,inVars,outVarsP);
        dG4p = subs(dG4,inVars,outVarsP);
        dG5p = subs(dG5,inVars,outVarsP);
        
        %%%% AMB - entire system @ New stance foot
        eqnHs0m = ...   %Before collision
            cross2d(G1m-P5m,m1*dG1m) + dq1m*I1 + ...
            cross2d(G2m-P5m,m2*dG2m) + dq2m*I2 + ...
            cross2d(G3m-P5m,m3*dG3m) + dq3m*I3 + ...
            cross2d(G4m-P5m,m4*dG4m) + dq4m*I4 + ...
            cross2d(G5m-P5m,m5*dG5m) + dq5m*I5;
        eqnHs0 = ...   %After collision
            cross2d(G1p-P0p,m1*dG1p) + dq1p*I1 + ...
            cross2d(G2p-P0p,m2*dG2p) + dq2p*I2 + ...
            cross2d(G3p-P0p,m3*dG3p) + dq3p*I3 + ...
            cross2d(G4p-P0p,m4*dG4p) + dq4p*I4 + ...
            cross2d(G5p-P0p,m5*dG5p) + dq5p*I5;
        
        
        %%%% AMB - new swing leg, torso, stance femer  @  stance knee
        eqnHs1m = ...   %Before collision
            cross2d(G1m-P4m,m1*dG1m) + dq1m*I1 + ...
            cross2d(G2m-P4m,m2*dG2m) + dq2m*I2 + ...
            cross2d(G3m-P4m,m3*dG3m) + dq3m*I3 + ...
            cross2d(G4m-P4m,m4*dG4m) + dq4m*I4;
        eqnHs1 = ...   %After collision
            cross2d(G2p-P1p,m2*dG2p) + dq2p*I2 + ...
            cross2d(G3p-P1p,m3*dG3p) + dq3p*I3 + ...
            cross2d(G4p-P1p,m4*dG4p) + dq4p*I4 + ...
            cross2d(G5p-P1p,m5*dG5p) + dq5p*I5;
        
        
        %%%% AMB - swing leg, torso  @ new hip
        eqnHs2m = ...   %Before collision
            cross2d(G3m-P2m,m3*dG3m) + dq3m*I3 + ...
            cross2d(G2m-P2m,m2*dG2m) + dq2m*I2 + ...
            cross2d(G1m-P2m,m1*dG1m) + dq1m*I1;
        eqnHs2 = ...   %After collision
            cross2d(G3p-P2p,m3*dG3p) + dq3p*I3 + ...
            cross2d(G4p-P2p,m4*dG4p) + dq4p*I4 + ...
            cross2d(G5p-P2p,m5*dG5p) + dq5p*I5;
        
        
        %%%% AMB - swing leg @ new hip
        eqnHs3m = ...   %Before collision
            cross2d(G1m-P2m,m1*dG1m) + dq1m*I1 + ...
            cross2d(G2m-P2m,m2*dG2m) + dq2m*I2;
        eqnHs3 = ...   %After collision
            cross2d(G4p-P2p,m4*dG4p) + dq4p*I4 + ...
            cross2d(G5p-P2p,m5*dG5p) + dq5p*I5;
        
        %%%% AMB - swing tibia @ new swing knee
        eqnHs4m = ...   %Before collision
            cross2d(G1m-P1m,m1*dG1m) + dq1m*I1;
        eqnHs4 = ...   %After collision
            cross2d(G5p-P4p,m5*dG5p) + dq5p*I5;
        
        
        %%%% Collect and solve equations:
        eqnHs = [...
            eqnHs0m - eqnHs0;
            eqnHs1m - eqnHs1;
            eqnHs2m - eqnHs2;
            eqnHs3m - eqnHs3;
            eqnHs4m - eqnHs4];
        [MM, FF] = equationsToMatrix(eqnHs,dqp);
        
        %%%% Compute gradients:
        tp = sym('tp','real');   %Initial trajectory time
        tm = sym('tm','real');   %Final trajectory time
        zBnd = [tp;qp;dqp;up;tm;qm;dqm;um];
        [m, mi, mz, mzi, mzd] = computeGradients(MM,zBnd,empty);
        [f, fi, fz, fzi, fzd] = computeGradients(FF,zBnd,empty);
        
        % Heel-strike
        matlabFunction(m, mi, f, fi,...   %dynamics
            mz, mzi, mzd, fz, fzi, fzd,...  %gradients
            'file','autoGen_cst_heelStrike.m',...
            'vars',{...
            'q1p','q2p','q3p','q4p','q5p',...
            'q1m','q2m','q3m','q4m','q5m',...
            'dq1m','dq2m','dq3m','dq4m','dq5m',...
            'm1','m2','m3','m4','m5',...
            'I1','I2','I3','I4','I5',...
            'l1','l2','l3','l4','l5',...
            'c1','c2','c3','c4','c5','empty'});
        
        % Collision velocity of the swing foot:
        cst = [-dP5p(2); dP5m(2)];  %Swing foot velocity before and after collision (negative sign is intentional, since output is constrained to be negative);
        cstJac = jacobian(cst,zBnd);  %Gradient
        matlabFunction(cst, cstJac,...
            'file','autoGen_cst_footVel.m',...
            'vars',{...
            'q1p','q2p','q4p','q5p',...
            'q1m','q2m','q4m','q5m',...
            'dq1p','dq2p','dq4p','dq5p',...
            'dq1m','dq2m','dq4m','dq5m',...
            'l1','l2','l4','l5'});
        
        % Step length and height constraint:
        stepLength = sym('stepLength','real');
        ceq = [P5m(1)-stepLength; P5m(2)];
        ceqJac = jacobian(ceq,zBnd);  %Gradient
        matlabFunction(ceq, ceqJac,...
            'file','autoGen_cst_steplength.m',...
            'vars',{...
            'q1m','q2m','q4m','q5m',...
            'l1','l2','l4','l5','stepLength'});
    end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Mechanical Energy                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    function mechanicalEnergy()
        disp('Deriving mechanical energy...')
        
        %%%% Energy:
        KineticEnergy = ...
            0.5*m1*dot(dG1,dG1) + 0.5*I1*dq1^2 + ...
            0.5*m2*dot(dG2,dG2) + 0.5*I2*dq2^2 + ...
            0.5*m3*dot(dG3,dG3) + 0.5*I3*dq3^2 + ...
            0.5*m4*dot(dG4,dG4) + 0.5*I4*dq4^2 + ...
            0.5*m5*dot(dG5,dG5) + 0.5*I5*dq5^2;
        PotentialEnergy = ...
            m1*g*G1(2) + ...
            m2*g*G2(2) + ...
            m3*g*G3(2) + ...
            m4*g*G4(2) + ...
            m5*g*G5(2);
        
        
        matlabFunction(KineticEnergy, PotentialEnergy,...
            'file','autoGen_energy.m',...
            'vars',{...
            'q1','q2','q3','q4','q5',...
            'dq1','dq2','dq3','dq4','dq5',...
            'm1','m2','m3','m4','m5',...
            'I1','I2','I3','I4','I5',...
            'l1','l2','l3','l4',...
            'c1','c2','c3','c4','c5',...
            'g'},...
            'outputs',{'KE','PE'});
        
    end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                          Contact Forces                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


    function contactForces()
        
        %%%% Contact Forces:
        eqnForce5 = w1 + w2 + w3 + w4 + w5 + Fx*i + Fy*j;
        eqnInertia5 = (m1+m2+m3+m4+m5)*ddG;
        [AA,bb] = equationsToMatrix(eqnForce5-eqnInertia5,[Fx;Fy]);
        ContactForces = AA\bb;
        
        matlabFunction(ContactForces(1),ContactForces(2),...
            'file','autoGen_contactForce.m',...
            'vars',{...
            'q1','q2','q3','q4','q5',...
            'dq1','dq2','dq3','dq4','dq5',...
            'ddq1','ddq2','ddq3','ddq4','ddq5',...
            'm1','m2','m3','m4','m5',...
            'l1','l2','l3','l4',...
            'c1','c2','c3','c4','c5',...
            'g'},...
            'outputs',{'Fx','Fy'});
        
        
    end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Write Kinematics Files                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    function kinematics()
        disp('Writing kinematics files...')
        
        
        P = [P1; P2; P3; P4; P5];
        Gvec = [G1; G2; G3; G4; G5];
        
        % Used for plotting and animation
        matlabFunction(P,Gvec,'file','autoGen_getPoints.m',...
            'vars',{...
            'q1','q2','q3','q4','q5',...
            'l1','l2','l3','l4','l5',...
            'c1','c2','c3','c4','c5'},...
            'outputs',{'P','Gvec'});
        
    end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Helper Functions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m, mi, mz, mzi, dim] = computeGradients(M,z,empty)
%
% This function computes the gradients of a matrix M with respect the the
% variables in z, and then returns both the matrix and its gradient as
% column vectors of their non-zero elements, along with the linear indicies
% to unpack them. It also simplifies m and mz.
%
% INPUTS:
%   M = [na, nb] = symbolic matrix
%   z = [nc, 1] = symbolic vector
%
% OUTPUTS:
%   m = [nd, 1] = symbolic vector of non-zero elements in M
%   mi = [nd, 1] = linear indicies to map m --> [na,nb] matrix
%   mz = [ne, 1] = symbolic vector of non-zero elements in Mz
%   mzi = [ne, 1] = linear indicies to map mz --> [na,nb,nc] array
%   dim = [3,1] = [na,nb,nc] = dimensions of 3d version of mz
%

[na, nb] = size(M);
nc = size(z,1);
M = simplify(M);

mz2 = jacobian(M(:),z);  %Compute jacobian of M, by first reshaping M to be a column vector
mz3 = reshape(mz2,na,nb,nc); %Expand back out to a three-dimensional array
mz3 = simplify(mz3);

% Extract non-zero elements to a column vector:
mi = find(M);
m = M(mi);
mzi = find(mz3);
mz = mz3(mzi); mz = mz(:);  %Collapse to a column vector
dim = [na,nb,nc];

% Pad any constant terms with "empty" to permit vectorization:
m = vectorizeHack(m, z, empty);
mz = vectorizeHack(mz, z, empty);

end



function x = vectorizeHack(x, z, empty)
%
% This function searches for any elements of x that are not dependent on
% any element of z. In this case, the automatically generated code will
% fail to vectorize properly. One solution is to add an array of zeros
% (empty) to the element.
%
% x = column vector of symbolic expressions
% z = column vector of symbolic variables
% z = symbolic variable, which the user will set equal to zero.
%

% Compute dependencies
g = jacobian(x,z);

% Check for rows of x with no dependence on z
[n,m] = size(g);
idxConst = true(n,1);
for i=1:n
    for j=1:m
        if ~isequal(sym(0),g(i,j))
            idxConst(i) = false;
            break;
        end
    end
end

% Add empty to those enteries
x(idxConst) = x(idxConst) + empty;

end





