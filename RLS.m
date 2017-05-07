function [P,l] = RLS(a,b,Q,R,K,gamma,lambda,numIter)
% This function estimates the parameters P using Recursive Least Square
% (RLS) method.
%
% Args:
%   a: System matrix (for simulation) (n x n)
%   b: Input matrix  (for simulation) (n x m)
%   Q: State weight matrix            (n x n)
%   R: Input weight matrix            (m x m)
%   K: Gain/Policy vector --> u = K*x (m x n)
%   gamma: Discount factor (for reinforcement learning)
%   lambda: Forgetting factor to the RLS [0.98 <---> 1] 
%   numIter: Number of iterations for the RLS algorithm
%
% Returns:
%   P: A vector of estimated parameters 
%   l: sqrt(size(P,1)), used to reshape P into a symmetric matrix 
% -------------------------------------------------------------------------


%TODO or REMOVE



end

