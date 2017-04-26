function [ P ] = calculateAnalyticalP(a,b,r,S)
% This function calculates the P matrix (Eq.58) analytically.
% 
% Args: 
%   a: System matrix
%   b: Input matrix
%   r: Prediction Horizon
%   S: Cost-to-go kernel matrix for [x ur]'
% 
% Returns: 
%   P: Cost-to-go kernel matrix for [x u]'
%-------------------------------------------------------------
    n = size(a,1);
    numInputs = size(b,2);
    Sxu = S(1:n, n+1:n+r*numInputs);
    Suu = S(n+1:n+r*numInputs,n+1:n+r*numInputs);
    G = -pinv(Suu)*Sxu';
    GL = G(1:numInputs,:);
    
    Gx = GL*a;
    for i=1:r-2
       Gx = vertcat(Gx, GL*((a+b*GL)^i)*a); 
    end
    
    Gu = GL*b;
    for i=1:r-2
       Gu = vertcat(Gu, GL*((a+b*GL)^i)*b); 
    end
    
    S11 = S(1:n,1:n);
    S12 = S(1:n,n+1:n+numInputs);
    S13 = S(1:n,n+numInputs+1:end);
    S21 = S(n+1:n+numInputs, 1:n);
    S22 = S(n+1:n+numInputs, n+1:n+numInputs);
    S23 = S(n+1:n+numInputs, n+numInputs+1:end);
    S31 = S(n+numInputs+1:end, 1:n);
    S32 = S(n+numInputs+1:end, n+1:n+numInputs);
    S33 = S(n+numInputs+1:end, n+numInputs+1:end);

    
    P =  [S11+Gx'*S31+S13*Gx+Gx'*S33*Gx,   S12+Gx'*S32+S13*Gu+Gx'*S33*Gu;
         S21+Gu'*S31+S23*Gx+Gu'*S33*Gx,   S22+Gu'*S32+S23*Gu+Gu'*S33*Gu]; 
end

