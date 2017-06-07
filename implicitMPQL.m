function [x,u_history] = implicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals,numIter, plotFlag)
% This function brings the system to zero using
% the implicit version of MPQL.
%
% Args:
%   a: system matrix
%   b: input  matrix
%   Q: state weigh matrix 
%   R: input weight matrix
%   r: predictive horizon 
%   gamma: discount factor
%   x_init: initial state of the system
%   inputVals: matrix of allowed discrete input values (each row representing an input)
%   numIter: number of iterations for the simulation
%   plotFlag: Whether to plot the result
%
% Returns:
%   x: The state values throughout the simulation
%   u_history: Control Signals throughout the simulation
%--------------------------------------------------------------------------

TEST_MODE = false;

    %shapes
    [n,m] = size(b);
    
    %---S matrix (choose one)---
    S = calculateAnalyticalS(a,b,r,gamma,Q,R);
%     S = calculateNumericalS(a,b,r,gamma,Q,R);
    %----------------------------------------
    
    %---P matrix (choose one)---
%     P_anas = calculateAnalyticalPs(a,b,r,S);
    P_ana = calculateAnalyticalP(a,b,r,S);
%     P_num_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
%     P_num_batch = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
    P = P_ana;
    %----------------------------------------
    
    % simulation
    numInputVals = size(inputVals,2);
    % close_loop
    x = zeros(n,numIter);
    x(:,1) = x_init; 
    Qs = zeros(1,numInputVals^m); %cost to go from taking each action
    u_history = zeros(m, numIter);
    [actionSpace, actionSpaceLength] = buildActionSpace(inputVals);
    for k=1:numIter
        for i=1:actionSpaceLength %action space
            xu = vertcat(x(:,k), actionSpace(:,i));
            Qs(i) =  xu'*P*xu;
        end
        [argvalue, argmin] = min(Qs);
        optimalInput = actionSpace(:,argmin); %choose the first value of the best input row values
        u_history(:,k) = optimalInput;
        x(:,k+1) = a*x(:,k) + b*optimalInput;
    end

    %plots
    if(plotFlag)
        % open_loop
        x_open = zeros(n, numIter);
        x_open(:,1) = x_init;
        for k=1:numIter
        x_open(:,k+1) = a*x_open(:,k) + b*zeros(m,1);
        end
        
        subplot(3,1,1)
        plot(0:numIter, x_open)
        title('Open-Loop')
        ylabel('state signals')
        xlabel('time step')

        subplot(3,1,2)
        plot(0:numIter, x)
        title(['Close-Loop with r=' num2str(r) ' gamma=' num2str(gamma) ''])
        ylabel('state signals')
        xlabel('time step')

        subplot(3,1,3)
        plot(0:numIter-1, u_history)
        title('Control Signal') 
        ylabel('control input')
        xlabel('time step')
    end
    
    
    
    
    if TEST_MODE
        %--- GP numerical (RLS or Batch) ---
        Pxu = P(1:n, n+1:end);
        Puu = P(n+1:end, n+1:end);
        GP_num = -inv(Puu)*Pxu';
        eigGP_num = eig(a+b*GP_num);
        %-------------------

        %--- GP from Policy Iteration ---
        Sxu = S(1:n, n+1:n+r*m);
        Suu = S(n+1:n+r*m,n+1:n+r*m);
        G = -pinv(Suu)*Sxu';
        GL = G(1:m,:);
        [P_PI,GP_PI_hist] = calculateOptimalP_PI(a,b,Q,R,r,gamma,GL,10);

        Pxu = P_PI(1:n, n+1:end);
        Puu = P_PI(n+1:end, n+1:end);
        GP_PI = -inv(Puu)*Pxu';
        eigGP_PI = eig(a+b*GP_PI);
        %------------------------------

        %--- GLQR ---
        G_LQR = -dlqr(a,b,Q,R);
        eigG_LQR = eig(a+b*G_LQR);
        %------------

        %--- Display Test Results---
        disp('ImplicitMPQL: Test Mode Results:')
        disp('Gains:')
        G_LQR
        GP_num
        GP_PI

        disp('Eigen Values:')
        eigG_LQR 
        eigGP_num
        eigGP_PI

        figure
        plot(GP_PI_hist)
        hold on
        plot(repmat(GP_num,[size(GP_PI_hist,1),1]))
        hold off
        title('GP\_num vs GP\_IP')
        xlabel('Iteration Steps (for Policy Iteration)')
        ylabel('Gain values')
        
        figure
        plot(GP_PI_hist)
        hold on
        plot(repmat(G_LQR,[size(GP_PI_hist,1),1]))
        hold off
        title('G\_LQR vs GP\_IP')
        xlabel('Iteration Steps (for Policy Iteration)')
        ylabel('Gain values')
    end
end

