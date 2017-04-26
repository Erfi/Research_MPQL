% Erfan Azad
% Date: 10 March 2017

clear all;

randn('seed', 1)
%Continuous

%======System 1========
% ac = [0 1; -10 0.1];
% bc = [0 1]';
%======================

%=====System 2=========
%two frequency system
ac = [0 0 1 0; 0 0 0 1;
     -8 2 -.5 4; 1 -1.25 2 -2.5];
bc = [0 0 1 0]';
%======================

%Discrete
dt = 0.05;
[a,b] = c2d(ac, bc, dt);
c = eye(2,2);
d = zeros(2,1);

%shapes
[n,~] = size(a);

%weight matrices
Q = 2*eye(4);
R = 1;

r = 3;          %predictive horizon
gamma = 0.5;    %discount factor

% analytical S
Strue = calculateNumericalS(a,b,r,gamma,Q,R);

% simulation
maxNumIter = 200;
inputVals = [12,-12]; %possible input values
x_init = [5,2,-5,2]'; %initial x

% close_loop
x = zeros(n,maxNumIter);
x(:,1) = x_init;  
[us, numInputRows] = getAllInputCombinations(inputVals, r); %matrix of all possible inputs for 'r' steps
vs = zeros(1,numInputRows); %array of cost to go values roccesponding to each row of input forn us 
u_history = zeros(1, maxNumIter);
for k=1:maxNumIter
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
x_open = zeros(n, maxNumIter);
x_open(:,1) = x_init;
for k=1:maxNumIter
    x_open(:,k+1) = a*x_open(:,k) + b*0;
end

%plots
if(1)
    subplot(3,1,1)
    plot(0:maxNumIter, x_open)
    title('Open-Loop')
    ylabel('state signals')
    xlabel('time step')
    
    subplot(3,1,2)
    plot(0:maxNumIter, x)
    title(['Close-Loop with r=' num2str(r) ''])
    ylabel('state signals')
    xlabel('time step')

    subplot(3,1,3)
    plot(0:maxNumIter-1, u_history)
    title('Control Signal') 
    ylabel('control input')
    xlabel('time step')
end