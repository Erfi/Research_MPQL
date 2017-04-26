function [ AnalyticalS ] = calculateAnalyticalS( a,b,r,gamma,Q,R )
%This function calculates the analytical S matrix
% needed to aclculate V(k) = [x_k, ur_k]T * S * [x_k, ur_k]
%
%Args:
%   a: system matrix
%   b: input  matrix
%   r: predictive horizon 
%   gamma: discount factor
%   Q: state weigh matrix 
%   R: input weight matrix
    
    [n,~] = size(a);
    [numStates, numInputs] = size(b);

    capQ = zeros(n*r,n*r);
    capGammaQ = zeros(n*r,n*r);
    for i=1:r
       capQ(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = Q;
       capGammaQ(n*(i-1)+1:n*i,n*(i-1)+1:n*i)= eye(n,n) * (sqrt(gamma))^(i-1);
    end
    Qgamma = capGammaQ*capQ*capGammaQ;

    capR = zeros(numInputs*r, numInputs*r);
    capGammaR = zeros(numInputs*r, numInputs*r);  
    for i=1:r
        capR(numInputs*(i-1)+1:numInputs*i, numInputs*(i-1)+1:numInputs*i) = R;
        capGammaR(numInputs*(i-1)+1:numInputs*i, numInputs*(i-1)+1:numInputs*i) = eye(numInputs, numInputs) * (sqrt(gamma))^(i-1);
    end
    Rgamma = capGammaR*capR*capGammaR;

    % p1 and p2 matrices 
    P1 = zeros(n*r, n);
    for i=1:r
        P1(n*(i-1)+1:n*i,:) = a^i;
    end

    P2 = zeros(numStates*r, numInputs*r);
    for i=1:r %(for row=i)
        row = zeros(numStates, numInputs*r);
        for j=1:i %(for col=j)
            row(:,numInputs*(j-1)+1:numInputs*j) = a^(i-j)*b;
        end
       P2(numStates*(i-1)+1:numStates*i,:) = row;
    end
    
    AnalyticalS = [P1'*Qgamma*P1, P1'*Qgamma*P2;
                   P2'*Qgamma*P1, Rgamma+P2'*Qgamma*P2];

end

