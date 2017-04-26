function [ action ] = chooseAction( position, grid, epsilon, maximize)
%This function is responsible to choose action at each
%state (x,y) <=> (col, rows) position based on an epsilon greedy scheme
    [~, ~, numActions] = size(grid);
    if(maximize)
        [~, Index] = max(grid(position(2), position(1),:));
        borderVal = -inf;
    else
        [~, Index] = min(grid(position(2), position(1),:));
        borderVal = inf;
    end
    
    if rand > epsilon
        %take optimal action
        action = Index;
    else
        action = randi([1,numActions]);
        % in case of choosing an action that crashes through the border
        % we have to reselect.
        while(grid(position(2), position(1), action) == borderVal)
            action = randi([1,numActions]);
        end
    end
end

