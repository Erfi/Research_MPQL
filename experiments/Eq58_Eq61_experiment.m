% This is an experiment to see if Eq.58 [calculateAnalyticalP()] will be
% the same as Eq.61/66 [CalculateNumericalP()] if we add the extra term 
% gamma*U(k+r) term to Eq.61/66.
% We are expecting to get the same answer when r is small and gamma is
% close to 1.

%------System (marginally stable)------
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 30;
    gamma=0.8;
%------------------
S = calculateNumericalS(a,b,r,gamma,Q,R);
P_analytical = calculateAnalyticalP(a,b,r,S);

%---modifying & calculating Eq.61/66---
[n,m] = size(b);
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL = G(1:m,:);
numIter = 2*(n+m)^2; %2 * number of equations nessessary (so we will have enough rank)

for i=1:numIter
    x = zeros(n,r+2); % need x at k, k+1, and k+r+1 so we will fill this array for each iteration
    u = zeros(m, r); % need u at k, k+r so we will fill this array for each iteration 
    x(:,1) = randn(n,1); % randomly initialize the state
    u(:,1) = GL*x(:,1) + randn; % randomly initialize the input signal
    for j=1:r+2-1
       x(:,j+1) = a*x(:,j) + b*u(:,j);
       u(:,j+1) = GL*x(:,j+1);
    end
    xu_k = vertcat(x(:,1), u(:,1));
    xu_kp1 = vertcat(x(:,2),u(:,2));
    U_k = x(:,2)'*Q*x(:,2) + u(:,1)'*R*u(:,1);
    U_kpr = x(:,r+2)'*Q*x(:,r+2) + u(:,r+1)'*R*u(:,r+1);
    
    LHS(i,:) = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
    RHS(i,:) = U_k - (gamma^r)*U_kpr;
end

P_stacked = pinv(LHS)*RHS;
m = sqrt(size(P_stacked,1));
P_numerical = reshape(P_stacked, [m,m]);
%--------------------------------------

r
gamma
P_analytical
P_numerical