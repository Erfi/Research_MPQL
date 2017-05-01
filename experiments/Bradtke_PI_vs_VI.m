% This script implements Policy Improvement and Value Iteration as
% described in 
% "Reinforcement Learning Applied to Linear Quadratic Regulation" by 
% [Steven J. Bradtke]
% and
% "Adaptive Linear Quadratic Control Using Policy Iteration" by 
% [Steven J. Bradtke et al]
%
% The claim is that Policy Iteration always converges to the optimal
% Q-function while Value Iteration does not necessarily converges to the
% optimal Q-function.
%----------------------------------------------------------------------

%------System (marginally stable)------
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 10;
    gamma=0.8; %LQR uses gamma = 1.0
%------------------

%TODO
