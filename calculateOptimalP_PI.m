function [ P, GP ] = calculateOptimalP_PI( a,b,Q,R,r,gamma,SorG_init,usingS,numIter)
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
%   SorG_init: initial stabilizing gain/policy or the S matrix
%   usingS: Flag indicating the previous argument is G_init or S
%   numIter: Number of policy iterations to perform
%
% Returns:
%   GP: Gains for each iteration. (m*numIter x n)
%--------------------------------------------------------------------------
    [n,m] = size(b);
    if usingS
        Sxu = SorG_init(1:n, n+1:n+r*m);
        Suu = SorG_init(n+1:n+r*m,n+1:n+r*m);
        G = -pinv(Suu)*Sxu';
        G_init = G(1:m,:); 
    else
        G_init = SorG_init;
    end
    
    GP = G_init;
    for i = 1:numIter
    %----Policy Evaluation----
        P = calculateNumericalP_RLS(a,b,Q,R,r,gamma,GP(end-(m-1):end,:),false);
%         P = calculateNumericalP(a,b,Q,R,r,gamma,GP(end-(m-1):end,:),false);
    %----Policy Improvement----
        P_xu = P(1:n, n+1:end);
        P_uu = P(n+1:end, n+1:end);
        newGain = -pinv(P_uu)*P_xu';
        GP = vertcat(GP,newGain);
    end
end

