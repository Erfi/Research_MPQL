function [ G,GL,GF ] = extractGainFromS(S,n,m)
% This function takes a matrix S (n+r*m x n+r*m)
% and returns gain and sungains G, GL and GF
%
% Args:
%   n: number of elements in state X -> (nx1)
%   m: number of input elements in u -> (mx1)
%---------------------------------------------
    
    Sxu = S(1:n,n+1:end);
    Suu = S(n+1:end,n+1:end);
    G = -inv(Suu)*Sxu';
    GL = G(1:m,:);
    GF = G(m+1:end,:);
end

