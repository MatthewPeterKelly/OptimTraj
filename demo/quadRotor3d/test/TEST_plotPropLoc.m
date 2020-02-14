disp('Running: TEST_plotPropLoc.m') 
figure 

%% Test 1 - zero position, along z axis
d_prop = 0.1 ; 
location = [0 0 0]' ; 
ax = [0 0 1]' ; 
plotPropLoc(d_prop, location, ax) ; 

%% Test 2 - zero position, along new axis
d_prop = 0.1 ; 
location = [0 0 0]' ; 
ax = [0 1 0]' ; 
plotPropLoc(d_prop, location, ax) ; 

%% Test 3 - non-zero position, along z axis
d_prop = 0.5 ; 
location = [0.5 0.5 0]' ; 
ax = [0 0 1]' ; 
plotPropLoc(d_prop, location, ax) ; 

%% Test 4 - non-zero position, along new axis
d_prop = 0.5 ; 
location = [0.5 0.5 1]' ; 
ax = [0 1 1]' ; 
plotPropLoc(d_prop, location, ax) ; 

%%
disp('TEST_plotPropLoc.m ran without error') 