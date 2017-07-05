% This is an experiment to see if Eq.58 [calculateAnalyticalP()] will be
% the same as Eq.61/66 [CalculateNumericalP()] if we add the extra term 
% gamma*U(k+r) term to Eq.61/66.
% We are expecting to get the same answer when r is small and gamma is
% close to 1.

clear all;
%--- System model ---
[a,b,C,D,Q,R,ac,bc] = getSystemModel(1);
%--------------------

[n,m] = size(b);
dt = 0.05;
r = 10;
gamma = 1;

S = calculateAnalyticalS(a,b,r,gamma,Q,R);
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL = G(1:m,:)                              % GL from S matrix



P_analytical = calculateAnalyticalP(a,b,r,S);
P_ana_xu = P_analytical(1:n, n+1:end);
P_ana_uu = P_analytical(n+1:end, n+1:end);
GP_ana = -inv(P_ana_uu)*P_ana_xu'          % GP_ana from analytical P Eq.58

%---OptimalP using PI---
% this should produce the GLQR (with enough iterations)
% note that if numIter (last argument) == 1 we get the 
% same result as the calculateNumericalP()
% or the expanded version below (without the extra term for RHS)

% [P_PI,GP_PI_hist] = calculateOptimalP_PI(a,b,Q,R,r,gamma,GL,10);
% P_PI_xu = P_PI(1:n, n+1:end);
% P_PI_uu = P_PI(n+1:end, n+1:end);
% GP_PI = -inv(P_PI_uu)*P_PI_xu'             % GP_PI from Optimal P using policy iteration
%-----------------------


%---modifying & calculating Eq.61/66---
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
conditionLHS = cond(LHS)
P_stacked = pinv(LHS)*RHS;
m = sqrt(size(P_stacked,1));
P_numerical = reshape(P_stacked, [m,m]);
%--------------------------------------

Pxu = P_numerical(1:n, n+1:end);
Puu = P_numerical(n+1:end, n+1:end);
GP_num = -inv(Puu)*Pxu'                     % GP_num from numerical P Eq.61/66

%-----LQR gain------
GLQR = -dlqr(a,b,Q,R)                       % GLQR from LQR (analytically)

%% plotting
numIter = 50;
% plot(GP_PI_hist)
% hold on
plot(repmat(GLQR(1,1),[numIter+1,1]))
hold on
plot(repmat(GP_num(1,1),[numIter+1,1]))
title('GLQR vs GP_PI vs GP_num')
hold off

%% Eigs
eigLQR = eig(a+b*GLQR)
eigGL = eig(a+b*GL)
eigGP_ana = eig(a+b*GP_ana)
eigGP_num = eig(a+b*GP_num)
eigGP_PI = eig(a+b*GP_PI)
