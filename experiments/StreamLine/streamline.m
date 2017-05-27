% This is an experiment to see if our algorithm S->P->repeat works for a
% situation that the data is being fed in real time.

clear all;
%------System (marginally stable)------
if(1)
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 10;
    gamma=1; %1.0 --> LQR uses gamma = 1.0
end
%------------------
%---INITIAL SETUP---
[n,m] = size(b);
numIter = 2000;
x_init = randn(n,1);
Gain = randn(m,n);
S = (1e-10)*randn(n+r*m, n+r*m);
P = (1e-10)*randn(n+m, n+m);
X = [x_init, zeros(n,numIter-1)];
U = zeros(m,numIter-1);
V_S = eye((n+r*m)^2)*1e8; %inverse corrolation matrix for S update
V_P = eye((n+m)^2)*1e8;   %inverse corrolation matrix for P update
%-------------------

for k = 1:numIter %Main Loop
    X_k = X(:,k);
    U_k = Gain*X_k + 0.1*randn(m,1); %---PLUS EXPLORATION---
    X_kp1 = simStep(a,b,U_k,X_k);
    %- Store state and inputs - 
    U(:,k) = U_k;
    X(:,k+1) = X_kp1;
    U(:,k+1) = Gain*X_kp1 + 0.1*randn(m,1); %---PLUS EXPROLATION---
    
    if(k > r) % needed r (+1) steps before we can continue with S update
        %----Update S-----
        % we are using datapoints that are shifted back by r
        Ur_k = reshape(U(:,k-r:k-1), [r*m,1]);  %reshape to stack it vertically
        Ur_kp1 = reshape(U(:,k-r+1:k), [r*m,1]);%reshape to stack it vertically
        XUr_k = vertcat(X_k, Ur_k);
        XUr_kp1 = vertcat(X_kp1, Ur_kp1);
        Util_k = X(:,k-r+1)'*Q*X(:,k-r+1) + U(:,k-r)'*R*U(:,k-r);
        Util_kpr = X(:,k+1)'*Q*X(:,k+1) + U(:,k)'*R*U(:,k);
        LHS = kron(XUr_k', XUr_k') - gamma*kron(XUr_kp1', XUr_kp1');
        RHS = Util_k - (gamma^(r))*Util_kpr;
        %- RLS -
        S_stacked = S(:);
        lambda = 1; % forgetting factor for the RLS
        V_S = (1/lambda)*(V_S - (((V_S*LHS')*(LHS*V_S))/(1+LHS*V_S*LHS')));
        S_stacked = S_stacked + V_S*LHS'*(RHS - LHS*S_stacked);
        %- update S-
        l = n+r*m;
        S = reshape(S_stacked,[l,l]);
        %-----------------

        %----Extract GL----
        [~,GL,~] = extractGainFromS(S,n,m);
        %------------------

        %----Update P------
        XU_k = vertcat(X(:,k), U(:,k));
        XU_kp1 = vertcat(X(:,k+1), GL*X(:,k+1));
        Util_k = X(:,k+1)'*Q*X(:,k+1) + U(:,k)'*R*U(:,k);
        LHS = kron(XU_k', XU_k') -gamma*kron(XU_kp1', XU_kp1');
        RHS = Util_k;
        %-RLS-
        P_stacked = P(:);
        lambda = 1; % forgetting factor for the RLS
        V_P = (1/lambda)*(V_P - (((V_P*LHS')*(LHS*V_P))/(1+LHS*V_P*LHS')));
        P_stacked = P_stacked + V_P*LHS'*(RHS - LHS*P_stacked);
        %-update P-
        l = n+m;
        P = reshape(P_stacked, [l,l]);
        %------------------

        %----Update Gain----
        GP = extractGainFromP(P,n);
        Gain = GP;
        %-------------------
    end
    Gain
end
%% -- PLOTS --

subplot(2,1,1);
plot(X');
title('X history')
xlabel('iteration number');
ylabel('X')
subplot(2,1,2);
plot(U)
title('Input history')
xlabel('iteration number');
ylabel('input u')
 

eig(a+b*Gain)
Glqr = -dlqr(a,b,Q,R)
eig(a+b*Glqr)
