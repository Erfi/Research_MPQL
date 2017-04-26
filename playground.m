clear all
randn('seed', 1)

% Us = getAllInputCombinations([-1,1],3);

%========System 1=========
if (0)
    ac = [0 1; -50 -0.3];
    bc = [0 1]';
    %Discrete
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    %weight matrices
    Q = 2*eye(2);
    R = 1;
    r = 3;          %predictive horizon
    gamma = 0.5;    %discount factor
    x_init = [2;-1];
    inputVals = [-1,1];
end
%=========================

%=====System 2 (walking robot)======
if(1)
   a=eye(2,2);
   b=eye(2,2);
%    dt=0.05;
%    [a,b] = c2d(a,b,dt);
   c=eye(2,2);
   d=zeros(2,2);
   Q=eye(2,2); 
   R=eye(2,2);
   gamma = 0.5;
   r=3;
   x_init=[20,10];
   inputVals = [-1,0,1;
                -1,0,1];
end
%===================================

%=======System3(two frequency)==========
if (0)
    ac = [0 0 1 0; 0 0 0 1;
          -8 2 -.5 4; 1 -1.25 2 -2.5];
    bc = [0 0 1 0]';
    %Discrete
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    %weight matrices
    Q = 2*eye(4);
    R = 1;
    r = 3;          %predictive horizon
    gamma = 0.5;    %discount factor
    inputVals = [12,-12]; %possible input values
    x_init = [5,2,-5,2]'; %initial x
end
%=======================================

%========TESTING=======
% S_true = calculateAnalyticalS(a,b,r, gamma, Q, R);
% S_estimate = calculateNumericalS(a,b,r,gamma,Q,R);
%---------------------------------------------------
% actionSpace = buildActionSpace(inputVals);
%---------------------------------------------------
% explicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals,300); %for now only works with single input (generalize getAllInputCombinations())
%---------------------------------------------------
[xImp, uImp] = implicitMPQL(a,b,Q,R,r,gamma,x_init,inputVals, 50);
figure;
scatter(xImp(1,:), xImp(2,:))
%---------------------------------------------------







