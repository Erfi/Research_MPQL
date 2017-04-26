% Erfan Azad
% Date: 10 March 2017

clear all;

randn('seed', 1)
%Continuous
% ac = [0 1; -10 0.1];
% bc = [0 1]';

%two frequency system
ac = [0 0 1 0; 0 0 0 1;
     -8 2 -.5 4; 1 -1.25 2 -2.5];
bc = [0 0 1 0]';

%Discrete
dt = 0.05;
[a,b] = c2d(ac, bc, dt);
c = eye(2,2);
d = zeros(2,1);

%shapes
[n,~] = size(a);
[numStates, numInputs] = size(b);

%weight matrices
Q = 2*eye(4);
R = 1;

r = 100;          %predictive horizon
gamma = 0.99;    %discount factor

% 
% maxNumIter = 400;
% inputVals = [13,8,3,0,-3,-8,-13]; %possible input values
% x_init = [5,2,-5,2]'; %initial x
% implicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals,maxNumIter,true);


% analytical S
Strue = calculateAnalyticalS(a,b,r,gamma,Q,R);
Sxu = Strue(1:n, n+1:n+r*numInputs);
Suu = Strue(n+1:n+r*numInputs,n+1:n+r*numInputs);
Gain = -pinv(Suu)*Sxu';
% simulation
maxNumIter = 400;
inputVals = [13,8,3,0,-3,-8,-13]; %possible input values
numInputVals = size(inputVals,2);
x_init = [5,2,-5,2]'; %initial x
% x_init = [5,2]';

% close_loop
x = zeros(n,maxNumIter);
x(:,1) = x_init; 
vs = zeros(1,numInputVals^numInputs); %cost to go from taking each action
u_history = zeros(numInputs, maxNumIter);
for k=1:maxNumIter
    us = Gain*x(:,k); % Optimal action sequence
    for i=1:numInputVals
        us(1:numInputs,1) = inputVals(1:numInputs,i); % change the first u(k) with one of the desrete possible actions
        xu = vertcat(x(:,k), us);
        vs(i) =  xu'*Strue*xu;
    end
    [argvalue, argmin] = min(vs);
    optimalInput = inputVals(argmin); %choose the first value of the best input row values
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
    title(['Close-Loop with r=' num2str(r) ' gamma=' num2str(gamma) ''])
    ylabel('state signals')
    xlabel('time step')

    subplot(3,1,3)
    plot(0:maxNumIter-1, u_history)
    title(['Control Signal -- inputs=[' num2str(inputVals) ']']) 
    ylabel('control input')
    xlabel('time step')
end