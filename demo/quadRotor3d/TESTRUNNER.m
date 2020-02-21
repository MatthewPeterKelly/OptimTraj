% runs all tests in 'test' folder

clc; clear;
addpath ./test ./utilities

tests = dir('./test/*TEST*.m') ;

disp('running TESTRUNNER') 
for i = 1:numel(tests)
    run(tests(i).name)
    clearvars -except tests
end

%% 
disp('TESTRUNNER ran without error') 