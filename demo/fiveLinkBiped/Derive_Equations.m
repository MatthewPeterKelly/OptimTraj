%%%% Derive Equations - Five Link Biped Model %%%%
%
% This script derives the equations of motion, as well as some other useful
% equations (kinematics, contact forces, ...) for the five-link biped
% model. 
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

%%%% Link state and derivatives
syms q1 q2 q3 q4 q5 'real' % Absolute link orientations
syms dq1 dq2 dq3 dq4 dq5 'real' % Absolute link angular rates
syms ddq1 ddq2 ddq3 ddq4 ddq5 'real' % Absolute link angular accelerations

%%%% System inputs (controls)
syms u1 'real' % Torque acting on stance leg tibia from ground
syms u2 'real' % torque acting on stance leg femur from stance leg tibia
syms u3 'real' % torque acting on torso from the stance leg femur
syms u4 'real' % torque acting on swing leg femur from torso
syms u5 'real' % torque acting on swing leg tibia from swing leg femur

%%%% Physical paramters 
syms m1 m2 m3 m4 m5 'real' % Link masses
syms c1 c2 c3 c4 c5 'real' % center of mass distance from joint
syms l1 l2 l3 l4 l5 'real' % link length
syms I1 I2 I3 I4 I4 'real' % link moment of inertia about center of mass
syms g 'real'  %Acceleration due to gravity

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Write Kinematics Files                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

P = [P1; P2; P3; P4; P5];
G = [G1; G2; G3; G4; G5];
matlabFunction(P,G,'file','autoGen_getPoints.m',...
    'vars',{...
    'q1','q2','q3','q4','q5',...
    'l1','l2','l3','l4','l5',...
    'c1','c2','c3','c4','c5'},...
    'outputs',{'P','G'});















