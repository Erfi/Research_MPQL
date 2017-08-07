function [ G_MPC ] = getMPCGain(A,B,Q,R,r,gamma)
% This function calculates the MPC gain using the original version.
% The first m rows of: -inv(R + P2'*Q*P2)*P2'*Q*P1
% 
% Args: 
%   A: System Dynamic Matrix
%   B: System Input Matrix
%   Q: State Weight Matrix
%   R: Input Weight Matrix
%   r: Predictive Horizon
% 
% Returns:
%   G_MPC: MPC gain
%--------------------------------------------------------------------------
[n,m] = size(B);
%--- make P1 ---
P1 = [];
for i=1:r
    P1 = vertcat(P1, A^i);
end
%--- make P2 ---
col = [];
for i=1:r
   col =  vertcat(col, (A^(i-1))*B);
end
P2 = [];
for i=1:r
   rowEnd = (r-(i-1))*n;
   currentCol = vertcat(zeros((i-1)*n,m), col(1:rowEnd,:));
   P2 = horzcat(P2, currentCol);
end

%--- calculate Qgamma and Rgamma ---
QCell = repmat({Q}, 1, r);
capQ = blkdiag(QCell{:});

RCell = repmat({R}, 1, r);
capR = blkdiag(RCell{:});

%--- Incorporating the discount factor ---
if (gamma ~= 1)
    capGammaQ = zeros(n*r,n*r);
    for i=1:r
       capGammaQ(n*(i-1)+1:n*i,n*(i-1)+1:n*i)= eye(n,n) * (sqrt(gamma))^(i-1);
    end
    Qgamma = capGammaQ*capQ*capGammaQ;

    capGammaR = zeros(m*r, m*r);  
    for i=1:r
        capGammaR(m*(i-1)+1:m*i, m*(i-1)+1:m*i) = eye(m, m) * (sqrt(gamma))^(i-1);
    end
    Rgamma = capGammaR*capR*capGammaR;
else
   Qgamma = capQ;
   Rgamma = capR;
end
    
%--- Calculating the gain
G = -pinv(Rgamma+P2'*Qgamma*P2)*P2'*Qgamma*P1;
G_MPC = G(1:m,:);
end

