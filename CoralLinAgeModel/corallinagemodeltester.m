classdef corallinagemodeltester < matlab.unittest.TestCase
    
    methods (Test)
        function testIfRuns(testCase)
            numPeriods = 25;
            ppp = 12;
            x = 0:1/ppp:numPeriods;
            y = data(ppp,numPeriods);
            corallinagemodel([x;y]');
        end
        
        function testNumYears(testCase)
            numPeriods = 25;
            ppp = 12;
            x = 0:1/ppp:numPeriods;
            y = data(ppp,numPeriods);
            ts = corallinagemodel([x;y]');
            difference = abs(ts(end,1) - numPeriods);
            testCase.verifyLessThan(difference,1);
        end
        
        function testNumYears_noise_white_05(testCase)
            passRate = testNumYears_noise_white(5,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
                
        function testNumYears_noise_white_25(testCase)
            passRate = testNumYears_noise_white(25,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
                
        function testNumYears_noise_white_50(testCase)
            passRate = testNumYears_noise_white(50,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
                
        function testNumYears_noise_white_75(testCase)
            passRate = testNumYears_noise_white(75,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
                
        function testNumYears_noise_white_85(testCase)
            passRate = testNumYears_noise_white(85,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
                
        function testNumYears_noise_white_95(testCase)
            passRate = testNumYears_noise_white(95,100,1);
            %sprintf('pass rate = %d', passRate);
            testCase.assertGreaterThan(passRate,95);
        end
        
    end
    
end 

%--------------------------------------------------------------
% helper functions
%--------------------------------------------------------------

function a = period
    a = 2*pi;
end
function a = data(p,n)
    x = 0:(period/p):(n*period);
    a = sin(x);
end

function [passRate] = testNumYears_noise_white(noise, n, tolerance)

    passRate = 0;
    noise = noise/100;

    for i = 1:n
        numPeriods = 25;
        ppp = 12;
        x = 0:1/ppp:numPeriods;
        y = data(ppp,numPeriods) + noise*randn(size(x));
        ts = corallinagemodel([x;y]', 'numyears', numPeriods);
        difference = abs(ts(end,1) - numPeriods);
        passRate = passRate + (difference < tolerance);
    end

    passRate = 100 * passRate/n;
end