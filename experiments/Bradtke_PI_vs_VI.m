% This script implements Policy Improvement and Value Iteration as
% described in 
% "Reinforcement Learning Applied to Linear Quadratic Regulation" by 
% [Steven J. Bradtke]
% and
% "Adaptive Linear Quadratic Control Using Policy Iteration" by 
% [Steven J. Bradtke et al]
%
% The claim is that Policy Iteration always converges to the optimal
% Q-function while Value Iteration does not necessarily converges to the
% optimal Q-function.
%----------------------------------------------------------------------
clear all;
% warning('off','all');
%------System (marginally stable)------
if(0)
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    gamma=1; %LQR uses gamma = 1.0
    K = [13,-2]; %initial gain/policy
end
%-------------------------------------
%-----System Phan's Experiment--------
if(1)
    ac=[0 1;-100 -1];
    bc=[0 1]';
    dt = 0.01;
    [a,b]=c2d(ac,bc,dt);
    c=eye(2,2);
    d=zeros(2,1);
    Q=[1 0;0 50];
    R=1;
    gamma = 1.0;
    K=[0 -1];
    xic=[1 -2]';
end
%-------------------------------------



%======POLICY ITERATION=======
if(1)
disp('Policy Iteration...');
numIters_Policy = 1500;
numIters_RLS = 100;
[n,m] = size(b);
Hu_stacked = zeros((n+m)^2,1); %Q-function parameterization matrix
Hu_stacked_temp = zeros((n+m)^2,1); % A backup of the Hu, to fall back on in case the new Hu produces unstable gain
for k = 1:numIters_Policy
    V_k = eye((n+m)^2)*10; % initialize autocorrolation matrix (large values --> low initial confidence)
    x_k = xic;
    for t = 1:numIters_RLS
        u_k = K*x_k + randn(m,1)*0.5; % u_k plus a random exploration component 
        x_kp1 = a*x_k + b*u_k; % simulating the system/plant/environment
        u_kp1 = K*x_kp1 + rand(m,1)*0.5; % u_kp1 plus a random exploration component
        %updating Hu using RLS
        xu_k = vertcat(x_k, u_k);
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_k'*Q*x_k + u_k'*R*u_k;
        LHS = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
        RHS = U_k;
        %RLS
        lambda = 1.0; % forgetting factor for the RLS
        V_kp1 = (1/lambda)*(V_k - (((V_k*LHS')*(LHS*V_k))/(1+LHS*V_k*LHS')));
        Hu_stacked_new = Hu_stacked_temp + V_kp1*LHS'*(RHS - LHS*Hu_stacked_temp);
        %update
        x_k = x_kp1;
        V_k = V_kp1;
        Hu_stacked_temp = Hu_stacked_new;
    end
    %Policy Improvement Step 
    l = sqrt(size(Hu_stacked_temp,1));
    Hu = reshape(Hu_stacked_temp, [l,l]);
    Hu_12 = Hu(1:n,n+1:end);
    Hu_22 = Hu(n+1:end, n+1:end);
    K_temp = -inv(Hu_22)*Hu_12';
    if max(abs(eig(a+b*K_temp))) > 1
        disp(['Unstable Gain..Skipping...PolicyIterNum: ', num2str(k)]);
        Hu_stacked_temp = Hu_stacked;
    else
        Hu_stacked = Hu_stacked_temp;
        K = K_temp;
    end
end
l = sqrt(size(Hu_stacked,1));
Hu_PI = reshape(Hu_stacked, [l,l]);
K_PI = K;
end
%==========END OF POLICY ITERATION========

%==============VALUE ITERATION============
if(0)
disp('Value Iteration...');
numIters_RLS = 50000;
[n,m] = size(b);
Hu_stacked = zeros((n+m)^2,1); %Q-function parameterization matrix
V_k = eye((n+m)^2)*10; % initialize autocorrolation matrix (large values --> low initial confidence)
x_k = [1,-5]';
% for t = 1:numIters_RLS
t = 0;
while (t < 10000)
    if mod(t,1000)==0
        V_k = eye((n+m)^2)*1;
    end
    u_k = K*x_k + randn(m,1)*0.5; % u_k plus a random exploration component 
    x_kp1 = a*x_k + b*u_k; % simulating the system/plant/environment
    u_kp1 = K*x_kp1 + randn(m,1)*0.5;
    %updating Hu using RLS
    xu_k = vertcat(x_k, u_k);
    xu_kp1 = vertcat(x_kp1,u_kp1);
    U_k = x_k'*Q*x_k + u_k'*R*u_k;
    LHS = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
    RHS = U_k;
    %RLS
    lambda =1; % forgetting factor for the RLS
    V_kp1 = (1/lambda)*(V_k - (((V_k*LHS')*(LHS*V_k))/(1+LHS*V_k*LHS')));
    Hu_stacked_new = Hu_stacked + V_kp1*LHS'*(RHS - LHS*Hu_stacked);
    %Policy Improvement Step 
    l = sqrt(size(Hu_stacked_new,1));
    Hu = reshape(Hu_stacked_new, [l,l]);
    Hu_12 = Hu(1:n,n+1:end);
    Hu_22 = Hu(n+1:end, n+1:end);
    K_temp = -inv(Hu_22)*Hu_12';
    %update
    if max(abs(eig(a+b*K_temp))) > 1
        disp(['Unstable Gain..Skipping...ValueIterNum: ', num2str(t)]);
    else
        Hu_stacked = Hu_stacked_new;
        K = K_temp;
        x_k = x_kp1;
        V_k = V_kp1;
        t = t+1;
    end
end
l = sqrt(size(Hu_stacked,1));
Hu_VI = reshape(Hu_stacked, [l,l]);
K_VI = K;
end
%==========END OF VALUE ITERATION=========

%-------K_LQR------
disp('--- LQR results ---');
K_LQR = -dlqr(a,b,Q,R);
K_LQR
eig_LQR = eig(a+b*K_LQR)

disp('--- Policy Iteration results ---');
Hu_PI
K_PI
eig_PI = eig(a+b*K_PI)

% disp('--- Value Iteration results ---');
% Hu_VI
% K_VI
% eig_VI = eig(a+b*K_VI)





