% This script is for comparing the result of 'calculateNumericalP', which
% uses batch identification and 'calculateNumericalP_RLS' which uses
% recursive identification (RLS).
% We are expecting the results to be the same.
%-------------------------------------------------------------------------
clear all;

%------System (marginally stable)------
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 10;
    gamma=0.8; %1.0 --> LQR uses gamma = 1.0
%------------------

S = calculateAnalyticalS(a,b,r,gamma,Q,R);
P_batch = calculateNumericalP(a,b,Q,R,r,gamma,S, true)
P_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S, true)