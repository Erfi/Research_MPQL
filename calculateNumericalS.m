function [ S_estimate ] = calculateNumericalS(a,b,r,gamma, Q,R)
%This funciton calculates the matrix S numerically/iteratively.
% S is needed to aclculate V(k) = [x_k, ur_k]T * S * [x_k, ur_k]
%
%Args:
%   a: system matrix
%   b: input  matrix
%   r: predictive horizon 
%   gamma: discount factor
%   Q: state weigh matrix 
%   R: input weight matrix

    [n,m] = size(b);
    numIter = 2*(n+r*m)^2; %twice more than neccessary so that we will get enough rank
%     LHS = zeros(numIter, ((n+r*m)+(n+r*m)^2)/2);
    LHS = zeros(numIter, (n+r*m)^2);
    RHS = zeros(numIter,1);
    for i=1:numIter 
        x = zeros(n, r+2);
        x(:,1) = 2*randn(n,1);
        us = randn(m,r+1);
        % calculating x(k), x(k+1) and x(k+r)
        for k=1:r+1
            x(:,k+1) = a*x(:,k) + b*us(:,k);
        end
        % building xu_k, xu_k+1
        xu_k = x(:,1);
        xu_kp1 = x(:,2);
        for k=1:r
            startIndex = n+m*(k-1)+1;
            endIndex = n+m*k;
            xu_k(startIndex:endIndex,1) = us(:,k);
            xu_kp1(startIndex:endIndex,1) = us(:,k+1);
            %or the version below wich is slightly slower but more readable
%             xu_k = vertcat(xu_k, us(:,k));
%             xu_kp1 = vertcat(xu_kp1, us(:,k+1));
        end
        % calculating U(k), the cost/reward of kth and k+rth step
        U_k = x(:,2)'*Q*x(:,2) + us(:,1)'*R*us(:,1);
        U_kpr = x(:,r+2)'*Q*x(:,r+2) + us(:,r+1)'*R*us(:,r+1);
        % calculating the LHS and RHS using (stack and) kronecker operators
        LHS(i,:) = kron(xu_k', xu_k') - gamma*kron(xu_kp1', xu_kp1');

%         LHS_full = kron(xu_k', xu_k') - gamma*kron(xu_kp1', xu_kp1');
%         LHS(i,:) = reduce_symmetric_vector(LHS_full);
        
        RHS(i,1) = U_k - (gamma^r)*U_kpr;
    end
    S_estimate_stacked = pinv(LHS)*RHS;
    m = sqrt(size(S_estimate_stacked,1));
    S_estimate = reshape(S_estimate_stacked, [m,m]);
%     S_estimate = expand_symmetric_vector(S_estimate_stacked, n+r*m);
end

