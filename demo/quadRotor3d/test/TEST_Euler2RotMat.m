disp('running TEST_Euler2RotMat.m') 

%% single vector
eul = zeros(1,3) ; 
[R] = Euler2RotMat(eul) ; 

%% wide array
eul = zeros(10,3) ; 
[R] = Euler2RotMat(eul) ; 

%%
disp('TEST_Euler3RotMat.m ran without error') 
