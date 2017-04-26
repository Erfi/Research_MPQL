function [ x,u_history ] = explicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals, numIter)
% This function brings the system to zero using
% the explicit version of MPQL.
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
%
% Returns:
%   x: The state values throughout the simulation
%   u_history: Control Signals throughout the simulation
%

    %shapes
    [n,~] = size(a);
    [numStates, numInputs] = size(b);
    % S matrix
    Strue = calculateNumericalS(a,b,r,gamma,Q,R);

    % close_loop
    x = zeros(n,numIter);
    x(:,1) = x_init;  
    [us, numInputRows] = getAllInputCombinations(inputVals, r); %matrix of all possible inputs for 'r' steps
    vs = zeros(1,numInputRows); %array of cost to go values roccesponding to each row of input forn us 
    u_history = zeros(1, numIter);
    for k=1:numIter
        for i=1:numInputRows
            xu = vertcat(x(:,k),us(i,:)');
            vs(i) =  xu'*Strue*xu;
        end
        [argvalue, argmin] = min(vs);
        optimalInput = (us(argmin,1)); %choose the first value of the best input row values
        u_history(k) = optimalInput;
        x(:,k+1) = a*x(:,k) + b*optimalInput;
    end

    % open_loop
    x_open = zeros(n, numIter);
    x_open(:,1) = x_init;
    for k=1:numIter
        x_open(:,k+1) = a*x_open(:,k) + b*0;
    end

    %plots
    if(1)
        subplot(3,1,1)
        plot(0:numIter, x_open)
        title('Open-Loop')
        ylabel('state signals')
        xlabel('time step')

        subplot(3,1,2)
        plot(0:numIter, x)
        title(['Close-Loop with r=' num2str(r) ''])
        ylabel('state signals')
        xlabel('time step')

        subplot(3,1,3)
        plot(0:numIter-1, u_history)
        title('Control Signal') 
        ylabel('control input')
        xlabel('time step')
    end
end

