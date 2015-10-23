function p = getPhysicalParameters()
% p = getPhysicalParameters()
%
% This function returns the parameters struct for the five link biped, with
% parameters that correspond to the robot RABBIT, taken from the 2003 paper
% by Westervelt, Grizzle, and Koditschek: "Hybrid Zero Dynamics of Planar
% Biped Walkers"
%

%%%% MASS
tibia_mass = 3.2;  %kg;
femur_mass = 6.8;  %kg
torso_mass = 20;  %kg
p.m1 = tibia_mass;
p.m2 = femur_mass;
p.m3 = torso_mass;
p.m4 = femur_mass;
p.m5 = tibia_mass;

%%%% MOMENT OF INERTIA (about link CoM)
tibia_inertia = 0.93;  %kg-m^2
femur_inertia = 1.08;  %kg-m^2
torso_inertia = 2.22;  %kg-m^2
p.I1 = tibia_inertia;
p.I2 = femur_inertia;
p.I3 = torso_inertia;
p.I4 = femur_inertia;
p.I5 = tibia_inertia;

%%%% LENGTH
tibia_length = 0.4;  %m
femur_length = 0.4;  %m
torso_length = 0.625;  %m
p.l1 = tibia_length;
p.l2 = femur_length;
p.l3 = torso_length;
p.l4 = femur_length;
p.l5 = tibia_length;

%%%% DISTANCE (parent joint to CoM)
tibia_dist = 0.128;  %m
femur_dist = 0.163;  %m
torso_dist = 0.2;  %m
p.c1 = tibia_dist;
p.c2 = femur_dist;
p.c3 = torso_dist;
p.c4 = femur_dist;
p.c5 = tibia_dist;

%%%% MISC
p.g = 9.81;  %Gravity acceleration

end