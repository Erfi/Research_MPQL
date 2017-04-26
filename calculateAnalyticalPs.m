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
    n = size(a,1);
    numInputs = size(b,2);
    
    Sxu = S(1:n, n+1:n+r*numInputs);
    Suu = S(n+1:n+r*numInputs,n+1:n+r*numInputs);
    G = -pinv(Suu)*Sxu';
    GF = G(2:end,:);
    
    Szz = S(1:n+1,1:n+1);
    Szx = S(1:n+1, n+2:n+r*numInputs);
    Sxx = S(n+2:n+r*numInputs,n+2:n+r*numInputs);
    SG = [Szz, Szx*GF;
          GF'*Szx', GF'*Sxx*GF];

    SG11 = SG(1:n,1:n);
    SG12 = SG(1:n,n+1:n+numInputs);
    SG13 = SG(1:n,n+numInputs+1:end);
    SG21 = SG(n+1:n+numInputs, 1:n);
    SG22 = SG(n+1:n+numInputs, n+1:n+numInputs);
    SG23 = SG(n+1:n+numInputs, n+numInputs+1:end);
    SG31 = SG(n+numInputs+1:n+numInputs+n, 1:n);
    SG32 = SG(n+numInputs+1:n+numInputs+n, n+1:n+numInputs);
    SG33 = SG(n+numInputs+1:n+numInputs+n, n+numInputs+1:end);

    Ps = [SG11+SG13+SG31+SG33 SG12+SG32;
         SG21+SG23 SG22];
end

