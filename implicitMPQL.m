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
%

    %shapes
    [n,m] = size(b);
    % analytical S
    S = calculateAnalyticalS(a,b,r,gamma,Q,R);
    P = calculateAnalyticalP(a,b,r,S);
    % simulation
    numInputVals = size(inputVals,2);
    % close_loop
    x = zeros(n,numIter);
    x(:,1) = x_init; 
    vs = zeros(1,numInputVals^m); %cost to go from taking each action
    u_history = zeros(m, numIter);
    [actionSpace, actionSpaceLength] = buildActionSpace(inputVals);
    for k=1:numIter
        for i=1:actionSpaceLength %action space
            xu = vertcat(x(:,k), actionSpace(:,i));
            vs(i) =  xu'*P*xu;
        end
        [argvalue, argmin] = min(vs);
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
end

