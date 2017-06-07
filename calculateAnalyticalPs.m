function [ Ps ] = calculateAnalyticalPs(a,b,r,S)
% This funciton calculates the Ps matrix (Eq.47) analytically
% 
% Args: 
%   a: System matrix
%   b: Input matrix
%   r: Prediction Horizon
%   S: Cost-to-go kernel matrix for [x ur]'
% 
% Returns: 
%   Ps: Cost-to-go kernel matrix for [x u]'
%-------------------------------------------------------------
    [n,m] = size(b);
    
    
    Sxu = S(1:n, n+1:n+r*m);
    Suu = S(n+1:n+r*m,n+1:n+r*m);
    G = -pinv(Suu)*Sxu';
    GF = G(2:end,:);
    
    Szz = S(1:n+1,1:n+1);
    Szx = S(1:n+1, n+2:n+r*m);
    Sxx = S(n+2:n+r*m,n+2:n+r*m);
    SG = [Szz, Szx*GF;
          GF'*Szx', GF'*Sxx*GF];

    SG11 = SG(1:n,1:n);
    SG12 = SG(1:n,n+1:n+m);
    SG13 = SG(1:n,n+m+1:end);
    SG21 = SG(n+1:n+m, 1:n);
    SG22 = SG(n+1:n+m, n+1:n+m);
    SG23 = SG(n+1:n+m, n+m+1:end);
    SG31 = SG(n+m+1:n+m+n, 1:n);
    SG32 = SG(n+m+1:n+m+n, n+1:n+m);
    SG33 = SG(n+m+1:n+m+n, n+m+1:end);

    Ps = [SG11+SG13+SG31+SG33 SG12+SG32;
         SG21+SG23 SG22];
end

