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
    r = 100;
    gamma=1; %1.0 --> LQR uses gamma = 1.0
%------------------
[n,m] = size(b);
%------------------
S = calculateAnalyticalS(a,b,r,gamma,Q,R);
P_ana = calculateAnalyticalP(a,b,r,S);
P_batch = calculateNumericalP(a,b,Q,R,r,gamma,S, true);
P_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S, true);
%------------------
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL = G(1:m,:)
%------------------
P_ana_xu = P_ana(1:n, n+1:end);
P_ana_uu = P_ana(n+1:end, n+1:end);
G_Pa = -inv(P_ana_uu)*P_ana_xu'
%------------------
P_batch_xu = P_batch(1:n, n+1:end);
P_batch_uu = P_batch(n+1:end, n+1:end);
G_Pb = -inv(P_batch_uu)*P_batch_xu'
%------------------
P_RLS_xu = P_RLS(1:n, n+1:end);
P_RLS_uu = P_RLS(n+1:end, n+1:end);
G_Pr = -inv(P_RLS_uu)*P_RLS_xu'
%------------------
G_LQR = -dlqr(a,b,Q,R)
%------------------

%-----Analysis------
% From the result the following conclusion can be achieved:
% 1) Calculating the P matrix using batch form or using RLS results in the
% same P matrix (as expected)
% 
% 2) The resulting matrix P is the parameterization for the FIRST converged
% Q-function. Meaning that THIS IS NOT THE OPTIMAL Q-FUNCTION. This is why
% the gain value does not match with the LQR version. In order to achieve
% the optimal LQR gain the P matrix should be recalculated starting with
% the gain derived from itself in a Policy Iteration manner.
% -------------------
