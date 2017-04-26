% Erfan Azad <erfan@dartmouth.edu>
% Date: 6 April 2017
% This file contains simulation for a tabular Q-learning algorithm
%----------
dim = [20,20];
goal = [1,1];
maximizeFlag = false;
epsilon = 0.8;
gamma = 0.9;
alpha = 1.0;
numIter = 500;
plotFlag = true;
obstacles = [5,2;5,3;5,4;5,5;5,6;5,7;5,8;5,9;5,10;
             5,6,;6,6;7,6;8,6;
             8,7;8,8;8,9];

trainedGrid = trainQL(dim(1),dim(2),goal,maximizeFlag,obstacles, epsilon, gamma, alpha, numIter, plotFlag);
%TODO: make a function for doing a path finding given the trainedGrid and
%epsilon=0