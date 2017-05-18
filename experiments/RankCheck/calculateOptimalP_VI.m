function [ P,GP, x_hist, u_hist, LHS_hist, RHS_hist] = calculateOptimalP_VI(a,b,Q,R,r,gamma,G_init)
% This function [tries to] calculate the optimal P matrix using the Value
% Iteration scheme. NOTE: We are not expecting the correct answer here. 
% This function is written as an experiment for analysis. 
%--------------------------------------------------------------------------
    [n,m] = size(b);
    %---------------
    GP = G_init;
    P_stacked = zeros((n+m)^2,1); % initialize P (in stacked form) 
    V_k = eye((n+m)^2)*1e8; % initialize autocorrolation matrix (Large numbers suggest low confidance for the initial estimation)
    x_k = randn(n,1); % initialize x(k)
    diffAvg = 1000; % average difference in P_stacked for each iteration
    epsilon = 1e-3; % ending criteria
    counter = 1;
    while(diffAvg > epsilon)
        u_k = GP(counter,:)*x_k + randn(m,1); % u(k) with random/exploritory component
        xu_k = vertcat(x_k, u_k);
        x_kp1 = a*x_k + b*u_k; % this is a simulation by the enviroment (plant)
        u_kp1 = GP(counter,:)*x_kp1;%H*x_k; % The two are the same
        xu_kp1 = vertcat(x_kp1,u_kp1);
        U_k = x_kp1'*Q*x_kp1 + u_k'*R*u_k;
        LHS = kron(xu_k', xu_k') -gamma*kron(xu_kp1', xu_kp1');
        RHS = U_k;
        %-----history-----
        LHS_hist(counter,:) = LHS;
        RHS_hist(counter,:) = RHS;
        x_hist(counter,:) = x_k';
        u_hist(counter,:) = u_k';
        %-----------------
        %RLS
        lambda = 1; % forgetting factor for the RLS
        V_kp1 = (1/lambda)*(V_k - (((V_k*LHS')*(LHS*V_k))/(1+LHS*V_k*LHS')));
        P_stacked_new = P_stacked + V_kp1*LHS'*(RHS - LHS*P_stacked);
        %ending criteria
        diffAvg = 0.9*diffAvg + 0.1*(norm(P_stacked - P_stacked_new));
        %update
        x_k = x_kp1;
        V_k = V_kp1;
        P_stacked_temp = P_stacked_new;
        %bound parameter of P so they won't blow up
%         P_stacked = min(P_stacked, 1e305); %upper bound limit
%         P_stacked = max(P_stacked, -1e305); %lower bound
        %++++VALUE ITERATION++++
        l = sqrt(size(P_stacked_temp,1));
        PP = reshape(P_stacked_temp, [l,l]);
        PP_xu = PP(1:n,n+1:end);
        PP_uu = PP(n+1:end, n+1:end);
        GP(counter+1,:) = -inv(PP_uu)*PP_xu';
        if(max(abs(eig(a+b*GP(counter+1,:))))>1)
            GP(counter+1,:) = GP(counter,:);
            disp('skipping the unstable gain');
        else
            P_stacked = P_stacked_temp;
        end
        %+++++++++++++++++++++++ 
        counter = counter + 1;
    end
    counter
    l = sqrt(size(P_stacked,1));
    P = reshape(P_stacked, [l,l]);
    rankVI = rank(LHS_hist(1:end-2,:))
end

