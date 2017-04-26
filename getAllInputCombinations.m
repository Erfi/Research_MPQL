function [ Us, numInputCombos ] = getAllInputCombinations( inputVals, numSteps )
%This function returns an array of arrays,
%Us of size (size(inputVals,2)^numSteps x numSteps)
%containing all the input combinations of inputs
% e.g. inputVals = [0,1], numSteps = 2
% Us = [0 0;
%       0 1;
%       1 0;
%       1 1]
%
% Args:
%   inputVals: values that the input can take at each time step. e.g. [0,1]
%   numSteps: How long shoud each array of inputs be e.g. 3 --> [1,0,1]
    numInputVals = size(inputVals,2);
    numInputCombos = numInputVals^numSteps;
    Us = zeros(numInputCombos, numSteps);
    for i=0:numInputCombos-1
        U = dec2base(i, numInputVals, numSteps);
        for j=1:numSteps
            index = str2double(U(1,j)) + 1;
            Us(i+1,j) = inputVals(index);    
        end
    end
end

