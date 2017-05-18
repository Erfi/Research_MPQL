function [ GP ] = extractGainFromP(P,n)
% This function takes a matrix P (n+m x n+m)
% and returns gain GP.
%
% Args:
%   n: number of elements in state X -> (nx1)
%---------------------------------------------
    Pxu = P(1:n,n+1:end);
    Puu = P(n+1:end,n+1:end);
    GP = -inv(Puu)*Pxu';
end




