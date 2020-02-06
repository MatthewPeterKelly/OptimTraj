disp('running TEST_Euler2RotMat.m') 

%% single vector
eul = zeros(1,3) ; 
[R] = Euler2RotMat(eul)

%% single vectors
eul = zeros(10,3) ; 
[R] = Euler2RotMat(eul) ; 