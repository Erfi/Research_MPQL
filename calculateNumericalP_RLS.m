function [ P ] = calculateNumericalP_RLS(a,b,Q,R,r,gamma,SorGL,usingS)
% This funciton calculates the P matrix (Eq.66) numerically
% Using RLS (not batch) identification.
%
% Args:
%   a: System matrix (used for simulation)
%   b: Input Matrix  (used for simulation)
%   Q: State weight matrix
%   R: Input weight matrix
%   r: Prediction Horizon
%   gamma: Discount factor
%   SorGL: Cost-to-go kernel matrix for [x ur]' or the Gain extracted from
%       it
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
    P_stacked = zeros((n+m)^2,1); % initialize P (in stacked form) 
    V_k = eye((n+m)^2)*100; % initialize autocorrolation matrix (Large numbers suggest low confidance for the initial estimation)
    x_k = randn(n,1); % initialize x(k)
    diffAvg = 100; % average difference in P_stacked for each iteration
    epsilon = 1e-6; % ending criteria
    while(diffAvg > epsilon)
        u_k = GL*x_k + 0.01*randn; % u(k) with random/exploritory component 
        xu_k = vertcat(x_k, u_k);
        x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
        u_kp1 = GL*x_kp1;%H*x_k; % The two are the same
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
        LHS = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
        RHS = U_k;
        %RLS
        lambda = 0.99; % forgetting factor for the RLS
        V_kp1 = (1/lambda)*(V_k - (((V_k*LHS')*(LHS*V_k))/(1+LHS*V_k*LHS')));
        P_stacked_new = P_stacked + V_kp1*LHS'*(RHS - LHS*P_stacked);
        %ending criteria
        diffAvg = 0.9*diffAvg + 0.1*sum(sum(abs(P_stacked - P_stacked_new)));
        %update
        x_k = x_kp1;
        V_k = V_kp1;
        P_stacked = P_stacked_new; 
    end
    l = sqrt(size(P_stacked,1));
    P = reshape(P_stacked, [l,l]);
end

