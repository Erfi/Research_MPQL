function [ P ] = calculateNumericalP(a,b,Q,R,r,gamma,SorGL,usingS)
% This funciton calculates the P matrix (Eq.66) numerically
% Using batch (not recursive e.g. RLS) identification.
%
% Args:
%   a: System matrix (used for simulation)
%   b: Input Matrix  (used for simulation)
%   Q: State weight matrix
%   R: Input weight matrix
%   r: Prediction Horizon
%   gamma: Discount factor
%   SorGL: Cost-to-go kernel matrix for [x ur]' or the gain, GL, extracted
%       from it.
%   usingS: Flag to indicate whether the second last input is matrix S or
%       the gain that is extracted from it, GL. This is used to calculate P 
%       iteratively in a policy iteration manner with a new GL everytime. 
%
% Returns:
%   P: Cost-to-go kernel matrix for [x u]'
%---------------------------------------------------------
    [n,m] = size(b);
    if usingS
        Sxu = SorGL(1:n, n+1:n+r*m);
        Suu = SorGL(n+1:n+r*m,n+1:n+r*m);
        G = -pinv(Suu)*Sxu';
        GL = G(1:m,:); 
    else
        GL = SorGL;
    end
    
    numIter = 2*((n+m)^2); % 2 * number of equations nessessary (so we will have enough rank)
    exploration_coeff = max(max(abs(GL))); % 100% of the max of GL
    
    for i=1:numIter
        x_k = randn(n,1); % initialize x(k)
        u_k = GL*x_k + exploration_coeff*randn(m,1);
        xu_k = vertcat(x_k, u_k);
        x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
        u_kp1 = GL*x_kp1;
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
        
        LHS(i,:) = kron(xu_k', xu_k') - gamma*kron(xu_kp1', xu_kp1');
%         LHS_full = kron(xu_k', xu_k') - gamma*kron(xu_kp1', xu_kp1');
%         LHS_reduced = reduce_symmetric_vector(LHS_full);
%         LHS(i,:) = LHS_reduced; 
        RHS(i,:) = U_k;
        %update state
%         x_k = x_kp1; %works better with independent data points 
    end
    %---- checking LHS matrix's condition ----
%     condition_Number = cond(LHS)
%     rank_Number = rank(LHS)
    %-----------------------------------------
    P_stacked = pinv(LHS)*RHS;
    m = sqrt(size(P_stacked,1));
    P = reshape(P_stacked, [m,m]);
%     P = expand_symmetric_vector(P_stacked,n+m);
end
