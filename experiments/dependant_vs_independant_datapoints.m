% This script is an experiment to test whether using depenent datapoints
% e.g. x(k)-->x(k+1)-->x(k+2)... gives the same result as using independent
% datapoints such as x(k)=randn, x(k+1)=randn,....
% We are expecting to get the same results.
%-------------------------------------------------------------------------

%------System (marginally stable)------
ac=[0 1; -50 0];
bc=[0 1]';
dt = 0.05;
[a,b] = c2d(ac, bc, dt);
Q=eye(2);
R=1;
r = 10;
gamma=0.8; %1.0 --> LQR uses gamma = 1.0
%------------------
%------Global Variables-----------
[n,m] = size(b);
S = calculateAnalyticalS(a,b,r,gamma,Q,R);
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL = G(1:m,:);
numIter = 2*(n+m)^2; %2 * number of equations nessessary (so we will have enough rank)
%---------------------------------
%-------Dependent Datapoints------
x_k = randn(n,1); % initialize the state
for i=1:numIter  
    u_k = GL*x_k + 0.01*randn;%randn;
    xu_k = vertcat(x_k, u_k);
    x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
    u_kp1 = GL*x_kp1;%H*x_k; % The two are the same
    xu_kp1 = vertcat(x_kp1,u_kp1);
    U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
    LHS(i,:) = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
    RHS(i,:) = U_k;
    %update
    x_k = x_kp1;
end
P_stacked = pinv(LHS)*RHS;
m = sqrt(size(P_stacked,1));
P_dependent = reshape(P_stacked, [m,m]);
%---------------------------------

%-----Independent Datapoints------
for i=1:numIter  
    u_k = GL*x_k + 0.01*randn;%randn;
    xu_k = vertcat(x_k, u_k);
    x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
    u_kp1 = GL*x_kp1;%H*x_k; % The two are the same
    xu_kp1 = vertcat(x_kp1,u_kp1);
    U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
    LHS(i,:) = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
    RHS(i,:) = U_k;
    %update
    x_k = x_kp1;
end
P_stacked = pinv(LHS)*RHS;
m = sqrt(size(P_stacked,1));
P_independent = reshape(P_stacked, [m,m]);
%---------------------------------

P_independent
P_dependent