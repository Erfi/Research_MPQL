function [ reward, nextPosition, done ] = step( position, action, goal )
% takes the given action in the given position
% and outputs the reward, and next position as well as if
% the goal has been reached.
    switch action
        case 1 %N
            nextPosition = [position(1),position(2)-1];
        case 2 %E
            nextPosition = [position(1)+1,position(2)];
        case 3 %S
            nextPosition = [position(1),position(2)+1];
        case 4 %W
            nextPosition = [position(1)-1,position(2)];
        otherwise
            disp('Invalid action value')
    end
    
    %using RL way
    
%     if(nextPosition == goal)
%         done = true;
%         reward = 1;
%     else
%         done = false;
%         reward = -1;
%     end

    %using CT way
    if(nextPosition ==  goal)
        done = true;
    else
        done = false;
    end
    
    Q = eye(2);
    R = 1;
    x_kp1 = abs(goal-nextPosition)';
    reward = x_kp1'*Q*x_kp1 + 1*R*1; 
end

