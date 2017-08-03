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
        G_init = extractGainFromS(SorG_init,n,m); 
    else
        G_init = SorG_init;
    end
    
    GP = G_init;
    for i = 1:numIter
        new_initial_gain = GP(end-(m-1):end,:);
    %----Policy Evaluation----
%         P = calculateNumericalP_RLS(a,b,Q,R,r,gamma,new_initial_gain,false);
        P = calculateNumericalP(a,b,Q,R,r,gamma,new_initial_gain,false);
    %----Policy Improvement----
        newGain = extractGainFromP(P,n);
        GP = vertcat(GP, newGain);
    end
end

