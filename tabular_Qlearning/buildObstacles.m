function [ obstacleGrid ] = buildObstacles( grid, obsCoords, maximize )
% This function takes a nx2 matrix containing coordinates [row, col]
% of the obstacles on the grid and puts -inf or inf depending on 
% if we are maximizing or minimizing the objective function.
%
% Args:
%   grid: The Q table containing the 4 Q values assisiated with each
%         action.
%   obsCoords: coordinates (row, col) of each obsticle in form of a nx2
%              matrix.
%   maximize: A flag indicating whether we are maximizing reward or
%             minimizing cost.
%
% Returns:
%   A grid containing obsticles meaning that the value of all actions that 
%   will result in crashing into onsticles are set to -inf or inf depending
%   on the 'maximize' parameter making them practically unselectable. 
%--------------------------------------------------------------------------

    if(maximize)
        value = -inf;
    else
        value = inf;
    end
    for i=1:size(obsCoords,1)
        row = obsCoords(i,1);
        col = obsCoords(i,2);
        grid(row+1,col,1) = value;
        grid(row,col-1,2) = value;
        grid(row-1, col,3) = value;
        grid(row, col+1,4) = value;
    end
    obstacleGrid = grid;
end

