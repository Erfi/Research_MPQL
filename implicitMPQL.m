function [x,u_history,GP,P] = implicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals,numIter,inputSelectionMode,plotFlag)
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
%   inputSelectionMode: <1: Selects the discrete input closest to the 
%                           continuous input according to GP> 
%                       <2: Uses P matrix to calculate Q values for each
%                           discrete inputs, and selects the one with
%                           minimum Q value>
%   plotFlag: Whether to plot the result
%
% Returns:
%   x: The state values throughout the simulation
%   u_history: Control Signals throughout the simulation
%--------------------------------------------------------------------------
    %shapes
    [n,m] = size(b);
    
    %---S matrix (choose one)---
    S = calculateAnalyticalS(a,b,r,gamma,Q,R);
%     S = calculateNumericalS(a,b,r,gamma,Q,R);
%     S = calculateNumericalS_RLS(a,b,r,gamma,Q,R);
    %----------------------------------------
    
    %---P matrix (choose one)---
%     P = calculateAnalyticalPs(a,b,r,S);
%     P = calculateAnalyticalP(a,b,r,S);
%     P = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
%     P = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
%     P = calculateOptimalP_PI(a,b,Q,R,r,gamma,S,true,5);
    %----------------------------------------
    
    
    %--- experiment to see if we can just get rid of the S matrix ---
%     randomG = randn(m,n)*1e-3;
    [~,randomG,~] = extractGainFromS(S,n,m);
%     mag_eig_initial_gain = abs(eig(a+b*randomG))
    P = calculateOptimalP_PI(a,b,Q,R,r,1.0,randomG,false,6);
    %----------------------------------------------------------------

    
    %----Calculating GP (for the output of this function)----
    GP = extractGainFromP(P,n)
    %---------------------------------------------------------
    
    if (inputSelectionMode == 1)
        %--- input selection based on GP ---
        numInputsVals = size(inputVals,2);
        x = zeros(n,numIter);
        x(:,1) = x_init;
        u_history = zeros(m, numIter);
        for k = 1:numIter
            uc = GP * x(:,k); %optimal continuous input
            u_k = zeros(m,1);
            [minVal, minIndex]=min(abs(inputVals - uc),[], 2);
            for i=1:m
               u_k(i) =  inputVals(i,minIndex(i));
            end
            u_history(:,k) = u_k;
            x(:,k+1) = a*x(:,k) + b*u_k;
        end
        %-------------------------------------

    elseif (inputSelectionMode == 2)
        %--- input selection based on Q values (P matrix) ---
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
            optimalInput = actionSpace(:,argmin); %choose the first value of the best input column values
            u_history(:,k) = optimalInput;
            x(:,k+1) = a*x(:,k) + b*optimalInput;
        end
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
end