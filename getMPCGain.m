function [ G_MPC ] = getMPCGain(A,B,Q,R,r)
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
%--- calculate G_MPC
QCell = repmat({Q}, 1, r);
Q_mat = blkdiag(QCell{:});

RCell = repmat({R}, 1, r);
R_mat = blkdiag(RCell{:});

G = -inv(R_mat+P2'*Q_mat*P2)*P2'*Q_mat*P1;
G_MPC = G(1:m,:);
end

