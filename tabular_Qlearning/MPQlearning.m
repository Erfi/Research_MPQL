% Erfan Azad <erfan@dartmouth.edu>
% Date: 10 April 2017
% This file contains the simulation for tabular MPQlearning.
% The problem is to create the Q-table using the MPQlearning
% algorithm for the walking robot with no obsticle.
%----------------------------------
clear all;
close all;

a=eye(2,2);
b=eye(2,2);
%    dt=0.05;
%    [a,b] = c2d(a,b,dt);
Q=eye(2,2); 
R=eye(2,2);
gamma = 1;
r=50;

        
gridDimX = 3;
gridDimY = 3;
numActions = 4; %[N,E,S,W]

S = calculateAnalyticalS(a,b,r,gamma,Q,R);
P = newPmatrix(a,b,r,S);
grid = makeGrid(gridDimX,gridDimY,numActions,false);
actions = [0,1,0,-1;
           -1,0,1,0]; %[N,E,S,W]

for row=1:gridDimY
    for col=1:gridDimX
        x = [col,row]'; % interpret as [x,y] <=> [col, row]
        for k = 1:4
            if(grid(row,col,k) ~= inf)
                u = actions(:,k);
                xu = vertcat(x,u);
                grid(row,col,k) = xu'*P*xu;
            end
        end
    end
end




gmax = gridMinVal(grid);
surf(gmax)
title('using MPQlearning');