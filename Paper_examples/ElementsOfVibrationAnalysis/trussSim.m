clear all;
close all;
%---------Book Example (Elements of Vibration Analysis)--------
m = 2; %kg
L = 2; %m
massMatrixCoeff = m*L/(6*sqrt(2));
massMatrix = [4*sqrt(2) 0 0 0 sqrt(2) 0 0;
              0 4+4*sqrt(2) 0 sqrt(2) 0 2 0;
              0 0 4+4*sqrt(2) 0 sqrt(2) 0 2;
              0 sqrt(2) 0 4+6*sqrt(2) 0 sqrt(2) 0;
              sqrt(2) 0 sqrt(2) 0 4+6*sqrt(2) 0 sqrt(2);
              0 2 0 sqrt(2) 0 4+2*sqrt(2) 0;
              0 0 2 0 sqrt(2) 0 4+2*sqrt(2)] * massMatrixCoeff;
E = 69e9; %Pa - Young's Modulus
A = pi*(0.005)^2;%m^2 Cross Sectional Area of beams
stiffnessMatrixCoeff = E*A/(2*sqrt(2)*L);
stiffnessMatrix = [2*sqrt(2) 0 0 0 0 0 0;
                   0 1+2*sqrt(2) -1 0 0 -1 1;
                   0 -1 1+2*sqrt(2) 0 -2*sqrt(2) 1 -1;
                   0 0 0 1+4*sqrt(2) -1 -2*sqrt(2) 0;
                   0 0 -2*sqrt(2) -1 1+2*sqrt(2) 0 0;
                   0 -1 1 -2*sqrt(2) 0 1+2*sqrt(2) -1;
                   0 1 -1 0 0 -1 1] * stiffnessMatrixCoeff;

% System dynamic (Continuous)
Bf = [0 0 0 0 0 0 1]';
Ac = [zeros(7,7) eye(7);
      -inv(massMatrix)*stiffnessMatrix zeros(7,7)];
Bc = [zeros(7,1);
      inv(massMatrix)*Bf];
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 1];
D = 0;
Q = eye(14);
R = 1e-6;

%------------Simulation Run Flags--------------
runContinuousLQR = true;
runDiscreteLQR = true;
runMPC = true;
runImplicitMPQL = true;

%------------Simulation Variables--------------
dt = 0.001;
Time = 0:dt:4;
U = zeros(size(Time));
X0 = [0 0 0 0 0 0 1 0 0 0 0 0 0 0]';

%--------------Continuous-Time-----------------
if(runContinuousLQR)
% Open-Loop Simulation (continuous)
sys_ol = ss(Ac,Bc,C,D);

[Y,~,~] = lsim(sys_ol, U, Time, X0);
figure(1);
subplot(2,1,1);
plot(Time,Y);
title('Open-Loop Simulation LQR (continuous-time)');

% Closed-Loop Simulation Using LQR (continuous)
K = -lqr(sys_ol, Q, R);        %LQR gain
Ac_cl = Ac+Bc*K;               %closed-Loop system dynamics matrix
sys_cl = ss(Ac_cl,Bc,C,D);

[Y,~,~] = lsim(sys_cl, U, Time, X0);
subplot(2,1,2)
plot(Time,Y);
title('Close-Loop Simulation LQR (continuous-time)');

end
%-----------------Discrete-Time-----------------
if(runDiscreteLQR)
[A,B] = c2d(Ac, Bc, dt);

%Open_Loop Simulation Using LQR (discrete)
[Y,~]=dlsim(A,B,C,D,U,X0);
figure(2);
subplot(3,1,1);
plot(Time,Y);
title('Open-Loop Simulation LQR (discrete-time)');

%Closed_Loop Simulation Using LQR (discrete)
Kd = -dlqr(A,B,Q,R);
[X_hist, U_hist] = simulate(A,B,Kd,X0,length(Time));
subplot(3,1,2);
plot(Time,X_hist(end,:));
title('Closed-Loop Simulation LQR (discrete-time)');

%Closed-Loop control signal
subplot(3,1,3);
plot(Time(1,1:end-1), U_hist);
title('Control Signal (input) LQR (discrete-time)');
end
%-------------------------MPC----------------------------
if(runMPC)
    [A,B] = c2d(Ac, Bc, dt);
    [n,numInputs] = size(B);
    r = 3;
    gamma = 0.9;
    S = calculateAnalyticalS(A,B,r,gamma,Q,R);
    Sxu = S(1:n, n+1:n+r*numInputs);
    Suu = S(n+1:n+r*numInputs,n+1:n+r*numInputs);
    G = -pinv(Suu)*Sxu';
    GL = G(1:numInputs,:);
    A_cl = A+B*GL;

    %Open_Loop Simulation Using LQR (discrete)
    [Y,~]=dlsim(A,B,C,D,U,X0);
    figure(3);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPC (discrete-time)');

    %Closed_Loop Simulation Using LQR (discrete)
    [X_hist, U_hist] = simulate(A,B,GL,X0,length(Time));
    subplot(3,1,2);
    plot(Time,X_hist(end,:))
    title(['Closed-Loop simulation MPC (discrete-time) with r=',num2str(r)]);
    
    %Closed-Loop control signal
    subplot(3,1,3);
    plot(Time(1,1:end-1), U_hist)
    title('Control Signal (input) MPC (discrete-time)');
end
%--------------------MPQL(implicit)----------------------
if(runImplicitMPQL)
    [A,B] = c2d(Ac, Bc, dt);
    r = 3;
    gamma = 0.9;
    input_vals = [-50000:10000:50000,-10000:1000:10000,-1000:100:1000,-100:10:100];
    numIter = length(Time);
    [X_hist, U_hist] = implicitMPQL(A,B,Q,R,r,gamma,X0,input_vals, numIter, false);
    
    %Open_Loop Simulation
    [Y,X]=dlsim(A,B,C,D,U,X0);
    figure(4);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPQL (discrete-time)'); 
    %Closed_Loop Simulation
    subplot(3,1,2);
    plot(Time,X_hist(14,1:end-1))
    title(['Closed-Loop simulation MPQL (discrete-time) with r=', num2str(r)]);
    %Control signal
    subplot(3,1,3);
    plot(Time, U_hist);
    title('Control Input');
end