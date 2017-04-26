function [ actionSpace, actionSpaceLength ] = buildActionSpace( inputVals )
% This function creates the action space for a given 
% matrix of possible inputs of size (numInputs, numAllowedValues)
% 
% Args:
%   inputVals: Matrix representing allowed input values for each
%              input. e.g. [-1,0,1; -2,0,2] for a two input system
%
% Returns:
%   Action Space: for the example above... (numInputs x numAllowedInputs^numInputs)
%    -1    -1    -1     0     0     0     1     1     1
%    -2     0     2    -2     0     2    -2     0     2
%--------------------------------------------------------------------------

    [numRows, numCols] = size(inputVals);
    actionSpaceLength = numCols^numRows;
    loopArray = ones(1, numRows);
    actionSpace = zeros(numRows, actionSpaceLength);
    for n=1:numCols^numRows
        %making a sequence according to the current loopArray
        seq = [];
        for i=1:numRows
            seq=[seq, inputVals(i,loopArray(i))];
        end
        %memorize the action sequence
        actionSpace(:,n) = seq';
        %incrementing the loopArray
        for k=numRows:-1:1
           if(loopArray(k) < numCols)
               loopArray(k) = loopArray(k) + 1;
               break;
           else
               loopArray(k) = 1;
           end
        end
    end
end

