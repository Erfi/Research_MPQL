function [ P ] = calculateNumericalP(a,b,Q,R,r,gamma,S)
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
%   S: Cost-to-go kernel matrix for [x ur]'
%
% Returns:
%   P: Cost-to-go kernel matrix for [x u]'
%---------------------------------------------------------
    [n,m] = size(b);
    Sxu = S(1:n, n+1:n+r*m);
    Suu = S(n+1:n+r*m,n+1:n+r*m);
    G = -pinv(Suu)*Sxu';
    GL = G(1:m,:);
    numIter = 2*(n+m)^2; %2 * number of equations nessessary (so we will have enough rank)
    for i=1:numIter
        x_k = randn(n,1);
        u_k = GL*x_k + 0.01*randn;%randn;
        xu_k = vertcat(x_k, u_k);
        x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
        u_kp1 = GL*x_kp1;%H*x_k; % The two are the same
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
        LHS(i,:) = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
        RHS(i,:) = U_k;
    end
    P_stacked = pinv(LHS)*RHS;
    m = sqrt(size(P_stacked,1));
    P = reshape(P_stacked, [m,m]);
end

