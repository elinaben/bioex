
classdef countATGInFramesTest < matlab.unittest.TestCase
    methods (Test)
        %define a test that gets and ORF and permutates it
        
        %for each 00
        
        function testSimple(testCase)
            % Define inputs and expected output
            count_ATGS = zeros(1,180);
            expected = count_ATGS;
            expected(90)=1;
            orf = "ATGGGAACATTAATGGTAGTA";
            utr5 = "GTGGGAACATTAATGGTAGTA";
            % Call function under test
            result = countATGInFrames(utr5, orf, length(utr5), max_length, count_ATGS)
            
            % Verify output
            testCase.verifyEqual(result, expected);
        end
    end
end
