function [ P ] = calculateNumericalP_RLS(a,b,Q,R,r,gamma,SorGL,usingS)
% This funciton calculates the P matrix (Eq.66) numerically
% Using RLS (not batch) identification. This is basically a policy
% evaluations step.
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
    [n,m] = size(b)
    if usingS
        Sxu = SorGL(1:n, n+1:n+r*m);
        Suu = SorGL(n+1:n+r*m,n+1:n+r*m);
        G = -pinv(Suu)*Sxu';
        GL = G(1:m,:); 
    else
         GL = SorGL;
    end
    P_stacked = zeros((n+m)^2,1); % initialize P (in stacked form) 
    V_k = eye((n+m)^2)*1e8; % initialize autocorrolation matrix (Large numbers suggest low confidance for the initial estimation)
    x_k = randn(n,1); % initialize x(k)
    for t=1:((n+m)^2)*5
        u_k = randn(m,1) + GL*x_k; % u(k) with random/exploritory component 
        xu_k = vertcat(x_k, u_k);
        x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
        u_kp1 = GL*x_kp1;%H*x_k; % The two are the same
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
        LHS = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
        RHS = U_k;
        %RLS step
        [V_k, P_stacked] = rls_one_step(V_k,P_stacked, LHS, RHS);
        %update state
        x_k = x_kp1;
    end
    l = sqrt(size(P_stacked,1));
    P = reshape(P_stacked, [l,l]);
end

