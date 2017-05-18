function [ P, GP ] = calculateOptimalP_PI( a,b,Q,R,r,gamma,G_init,numIter)
% This function calculates the optimal P using Policy Iteration. Policy
% Iteration consists of doing 1) policy evaluation 2) policy Improvement in
% a cyclic manner.
%
% Args:
%   a: System matrix (used for simulation)
%   b: Input Matrix  (used for simulation)
%   Q: State weight matrix
%   R: Input weight matrix
%   r: Prediction Horizon
%   gamma: Discount factor
%   G_init: initial stabilizing gain/policy
%   numIter: Number of policy iterations to perform
%
% Returns:
%   GP: Gains for each iteration. (m*numIter x n)
%--------------------------------------------------------------------------
    [n,m] = size(b);
    GP = G_init;
    for i = 1:numIter
        %----Policy Evaluation----
        P = calculateNumericalP_RLS(a,b,Q,R,r,gamma,GP(i,:),false);
        %----Policy Improvement----
        P_xu = P(1:n, n+1:end);
        P_uu = P(n+1:end, n+1:end);
        GP(i+1,:) = -inv(P_uu)*P_xu';
    end
end

