function [ grid ] = trainQL( dimX, dimY, goal, maximize, obstacles, epsilon, gamma, alpha, numIter, plotFlag )
%
% Args:
%   dimx: x dimention of the grid (number of columns).
%   dimY: y dimention of the grid (number of rowa).
%   goal: Terminal state (x,y) <=> (col, row).
%   obstacles: An nx2 matrix containing the coordinate (row, col) of
%              obstacles. Use '[]' (without quotations) for no obstacles.
%   epsilon: Exploration rate [0 to 1].
%   gamma: Discount factor for future cost/reward.
%   alpha: Learning rate [set it to 1 for deterministic systems].
%   numIter: Number of iterations for training.
%   plotFlag: A true/false flag to ndicate whether we want to plot the
%             resulting grid values.
% 
% Returns:
%   grid: Grid containing the Q values for each action.
%--------------------------------------------------------------------------

    gridDimX = dimX;
    gridDimY = dimY;
    numActions = 4; %[N,E,S,W]

%     goal = goal; %[x,y] = [col, row]
    maximizeFlag = maximize;
%     obstacles = [5,2;5,3;5,4;5,5;5,6;5,7;5,8;5,9;5,10
%                  5,6,;6,6;7,6;8,6;
%                  8,7;8,8;8,9];
    grid = makeGrid(gridDimX,gridDimY,numActions,obstacles,maximizeFlag);
    done = false; 
%     epsilon = 0.9;
%     gamma = 0.45;
%     alpha = 1; %deterministic env

    for k = 1:numIter
        iter = k
        position = [randi([1,gridDimX]), randi([1,gridDimY])];
        if k==200
            position = [15,15];
        end
        history = [];
        history = [history, position];
        done = false;
        while(~done)
           action = chooseAction(position, grid, epsilon, maximizeFlag);
           [reward, nextPosition, done] = step(position, action, goal);
           grid = updateQ(position, grid, action, reward, nextPosition, alpha, gamma, maximizeFlag);
           position = nextPosition;
           history = vertcat(history, position);
        end
        epsilon = max(0.05, epsilon-0.001);
    end

    if plotFlag
        if(maximizeFlag)
            Qgrid = gridMaxVal(grid);
        else
            Qgrid = gridMinVal(grid);
        end
        surf(Qgrid)
        title('using Q-learning with reward = utility (minimizing cost)');
    end

%     figure;
%     scatter(history(:,1), history(:,2))
%     title('last path and obsticles')
%     hold on;
%     scatter(obstacles(:,2), obstacles(:,1), 'filled')
end

