function [ updatedGrid ] = updateQ(position, grid, action, reward, nextPosition, alpha, gamma, maximize)
% updates the grid (Q-table) based on Q-learning formula
    row = position(2);
    col = position(1);
    nextRow = nextPosition(2);
    nextCol = nextPosition(1);
    if(maximize)
        grid(row, col, action) = grid(row, col, action) + alpha*(reward + gamma*max(grid(nextRow, nextCol,:)) - grid(row, col, action));
    else
        grid(row, col, action) = grid(row, col, action) + alpha*(reward + gamma*min(grid(nextRow, nextCol,:)) - grid(row, col, action));
    end
    updatedGrid = grid;
end

