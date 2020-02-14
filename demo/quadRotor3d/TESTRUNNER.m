% runs all tests in 'test' folder

tests = dir('./test/*TEST*.m') ;

for i = 1:numel(tests)
    run(tests(i).name)
    clearvars -except tests
end

%% 
disp('TESTRUNNER ran without error') 