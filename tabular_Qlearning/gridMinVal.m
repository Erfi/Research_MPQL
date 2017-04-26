function [ minGrid ] = gridMinVal( grid )
%creates a grid with max values at each position
    [yDim, xDim, numActions] = size(grid);
    minGrid = zeros(yDim, xDim);
    for row=1:yDim
        for col=1:xDim
            minGrid(row, col) = min(grid(row, col,:));
        end
    end
end
