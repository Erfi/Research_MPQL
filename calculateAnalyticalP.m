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
    [n,m] = size(b);
    Sxu = S(1:n, n+1:n+r*m);
    Suu = S(n+1:n+r*m,n+1:n+r*m);
    G = -pinv(Suu)*Sxu';
    GL = G(1:m,:);
    
    Gx = GL*a;
    for i=1:r-2
       Gx = vertcat(Gx, GL*((a+b*GL)^i)*a); 
    end
    
    Gu = GL*b;
    for i=1:r-2
       Gu = vertcat(Gu, GL*((a+b*GL)^i)*b); 
    end
    
    S11 = S(1:n,1:n);
    S12 = S(1:n,n+1:n+m);
    S13 = S(1:n,n+m+1:end);
    S21 = S(n+1:n+m, 1:n);
    S22 = S(n+1:n+m, n+1:n+m);
    S23 = S(n+1:n+m, n+m+1:end);
    S31 = S(n+m+1:end, 1:n);
    S32 = S(n+m+1:end, n+1:n+m);
    S33 = S(n+m+1:end, n+m+1:end);
    
    if r==1
        S13 = 0;
        S23 = 0;
        S33 = 0;
        S31 = 0;
        S32 = 0;
    end
    
    P =  [S11+Gx'*S31+S13*Gx+Gx'*S33*Gx,   S12+Gx'*S32+S13*Gu+Gx'*S33*Gu;
         S21+Gu'*S31+S23*Gx+Gu'*S33*Gx,   S22+Gu'*S32+S23*Gu+Gu'*S33*Gu]; 
end

