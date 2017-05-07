function [ S ] = calculateNumericalS_RLS( a,b,r,gamma, Q,R )
% This funciton calculates the matrix S numerically using 
% Recursive Least Square method (Not Batch method).
% S is needed to aclculate V(k) = [x_k, ur_k]T * S * [x_k, ur_k]
%
% Args:
%   a: system matrix (Used for simulation)
%   b: input  matrix (Used for simulation)
%   r: predictive horizon 
%   gamma: discount factor
%   Q: state weigh matrix 
%   R: input weight matrix
%
% Retuns:
%   S: The kernel matrix mapping [x ur]' to the cost-to-go V(x, ur) 
%--------------------------------------------------------------

    [n,m] = size(b);
    S_stacked = zeros((n+r*m)^2,1);
    V_k = eye((n+r*m)^2)*100; % initialize autocorrolation matrix (Large numbers suggest low confidance for the initial estimation)
    diffAvg = 100; % average difference in S_stacked for each iteration
    epsilon = 1e-8; % ending criteria
    while(diffAvg > epsilon) 
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
        LHS = kron(xu_k', xu_k') - gamma*kron(xu_kp1', xu_kp1');
        RHS = U_k - (gamma^r)*U_kpr;
        %RLS
        lambda = 0.99; % forgetting factor for the RLS
        V_kp1 = (1/lambda)*(V_k - (((V_k*LHS')*(LHS*V_k))/(1+LHS*V_k*LHS')));
        S_stacked_new = S_stacked + V_kp1*LHS'*(RHS - LHS*S_stacked);
        %ending criteria
        diffAvg = 0.9*diffAvg + 0.1*sum(sum(abs(S_stacked - S_stacked_new)));
        %update
        V_k = V_kp1;
        S_stacked = S_stacked_new;
    end
    l = sqrt(size(S_stacked,1));
    S = reshape(S_stacked, [l,l]);
end


