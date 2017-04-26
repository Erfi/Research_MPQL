function [ maxGrid ] = gridMaxVal( grid )
%creates a grid with max values at each position
    [yDim, xDim, numActions] = size(grid);
    maxGrid = zeros(yDim, xDim);
    for row=1:yDim
        for col=1:xDim
            maxGrid(row, col) = max(grid(row, col,:));
        end
    end
end

