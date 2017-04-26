function [ grid ] = makeGrid( gridDimX,gridDimY,numActions, obstacles, maximize)
%Makes the grid and puts -inf for inappropriate action vals
% e.g. The top row are not allowed to take action 1 (North)
% actions = [1,2,3,4] --> [N,E,S,W]

%     grid = zeros(gridDimY,gridDimX,numActions);
    grid = rand(gridDimY,gridDimX,numActions)*1e-5;
    
    if(maximize)
        grid(1,:,1) = -inf;
        grid(gridDimY,:,3) = -inf;
        grid(:,1,4) = -inf;
        grid(:,gridDimX,2) = -inf;

    else
        grid(1,:,1) = inf;
        grid(gridDimY,:,3) = inf;
        grid(:,1,4) = inf;
        grid(:,gridDimX,2) = inf;

    end
    
    grid = buildObstacles(grid, obstacles, maximize);
end

