% This is a sumulation using the matrices from the paper
% "Extracting Physical Parameters of Mechanical Models From 
% Identified State-Space Representation" by M. De Angelis et al.
%-----------------------------------------------------------------
clear all;
close all;

format short

[A,B,Cc,Dc,Q,R,Ac,Bc] = getSystemModel(4);
load('trussmodelforErfan.mat','dt');
% dt = 0.05;
[n,m] = size(B);

%--- Checking for controlibility ---
Co = ctrb(A,B);
if rank(Co) == n 
    disp('System is controlable')
else
    disp('System is not controllable')
end
%-----------------------------------
%--- Checking for observability  ---
Ob = obsv(A,Cc);
if rank(Ob) == n
    disp('System is observable')
else
    disp('System is not observable')
end
%-----------------------------------

%---- Use to decide on the sampling frequency dt ---
if(0)
%     figure(100)
%     subplot(2,1,1)
%     impulse(Ac,Bc, Cc, Dc) % Continuous Impulse Response
%     title('continuous')
%     subplot(2,1,2)
%     dimpulse(A,B,Cc,Dc);
%     title('discrete')
    maxFreqAc = max(imag(eig(Ac)))/(2*pi);
    idealSampleRate = 1.0/(6*maxFreqAc)
end
%----------------------------------------------------

%------------Simulation Run Flags--------------------
runContinuousLQR = false;
runDiscreteLQR =   true;
runMPC_Original =  true;
runMPC =           true;
runImplicitMPQL =  true;
%---------------------------------------------------

%------------Simulation Variables-------------------
Time = 0:dt:2;
U = zeros(size(Time,2),m); %single input
X0 = ones(1,n);
%----------------------------------------------------

%--------------Continuous-Time LQR-------------------
if(runContinuousLQR)
% Open-Loop Simulation (continuous)
sys_ol = ss(Ac,Bc,Cc,Dc);

[Yo,~,~] = lsim(sys_ol, U, Time, X0);
figure(1);
subplot(2,1,1);
plot(Time,Yo);
title('Open-Loop Simulation LQR (continuous-time)');

% Closed-Loop Simulation Using LQR (continuous)
K = -lqr(sys_ol, Q, R);        %LQR gain
Ac_cl = Ac+Bc*K;               %closed-Loop system dynamics matrix
sys_cl = ss(Ac_cl,Bc,Cc,Dc);

[Yc,~,~] = lsim(sys_cl, U, Time, X0);
subplot(2,1,2)
plot(Time,Yc);
title('Close-Loop Simulation LQR (continuous-time)');
end
%-----------------------------------------------------

%--------------Discrete-Time LQR----------------------
if(runDiscreteLQR)

%Open_Loop Simulation Using LQR (discrete)
[Y,~]=dlsim(A,B,Cc,Dc,U,X0);
figure(2);
subplot(3,1,1);
plot(Time,Y);
title('Open-Loop Simulation LQR (discrete-time)');

%Closed_Loop Simulation Using LQR (discrete)
Kd = -dlqr(A,B,Q,R)
[X_hist, U_hist] = simulate(A,B,Kd,X0,length(Time));
subplot(3,1,2);
plot(Time,Cc*X_hist);
title('Closed-Loop Simulation LQR (discrete-time)');

%Closed-Loop control signal
subplot(3,1,3);
plot(Time(1,1:end-1), U_hist);
title('Control Signal (input) LQR (discrete-time)');
end
%------------------------------------------------------

%------------------MPC (Original)----------------------
if (runMPC_Original)
    r = 50;
    G_MPC = getMPCGain(A,B,Q,R,r)

    %Open_Loop Simulation Using LQR (discrete)
    [Y,~]=dlsim(A,B,Cc,Dc,U,X0);
    figure(3);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPC\_Original');

    %Closed_Loop Simulation Using LQR (discrete)
    [X_hist, U_hist] = simulate(A,B,G_MPC,X0,length(Time));
    subplot(3,1,2);
    plot(Time,Cc*X_hist)
    title(['Closed-Loop simulation MPC\_Original with r=',num2str(r)]);
    
    %Closed-Loop control signal
    subplot(3,1,3);
    plot(Time(1,1:end-1), U_hist)
    title('Control Signal (input) MPC\_Original');
end

%------------------------------------------------------

%------------------MPC (from S)------------------------
if(runMPC)
    r = 50;
    gamma = 0.85;
    S = calculateAnalyticalS(A,B,r,gamma,Q,R);
    [~,GL,~] = extractGainFromS(S,n,m)

    %Open_Loop Simulation Using LQR (discrete)
    [Y,~]=dlsim(A,B,Cc,Dc,U,X0);
    figure(4);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPC (discrete-time)');

    %Closed_Loop Simulation Using LQR (discrete)
    [X_hist, U_hist] = simulate(A,B,GL,X0,length(Time));
    subplot(3,1,2);
    plot(Time,Cc*X_hist)
    title(['Closed-Loop simulation MPC (discrete-time) with r=',num2str(r),' gamma=', num2str(gamma)]);
    
    %Closed-Loop control signal
    subplot(3,1,3);
    plot(Time(1,1:end-1), U_hist)
    title('Control Signal (input) MPC (discrete-time)');
end
%--------------------------------------------------------

%--------------------MPQL(implicit)----------------------
if(runImplicitMPQL)
    
    r = 50;
    gamma = 0.86;
  
    input_vals = [-500000:100000:500000,-10000:1000:10000,-1000:100:1000,-100:10:100];
    input_vals = [-40000,40000,-30000,30000,-20000,20000,-10000,10000,-1000,1000,-100,100,0];
    input_vals = [-5000, 5000, -500,500, -50, 50, 0];
    input_vals = [-5,5,-4,4,-2,2,-1:0.2:1];
    input_vals = repmat(input_vals, [m,1]);
        
    numIter = length(Time);
    [X_hist_MPQL, U_hist_MPQL, GP] = implicitMPQL(A,B,Q,R,r,gamma,X0,input_vals, numIter,1, false);
    GP
   
    %Open_Loop Simulation
    [Y,X]=dlsim(A,B,Cc,Dc,U,X0);
    figure(5);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPQL (discrete-time)'); 
    %Closed_Loop Simulation
    subplot(3,1,2);
    plot(Time,Cc*X_hist_MPQL(:,1:end-1))
    title(['Closed-Loop simulation MPQL (discrete-time) with r=', num2str(r),' gamma=', num2str(gamma)]);
    %Control signal
    subplot(3,1,3);
    plot(Time, U_hist_MPQL);
    title('Control Input');
end
%---------------------------------------------------------

%% Looking at the magnitude of the eigenValues of GP:
format short
eigMag = abs(eig(A+B*GP))
% eigMagLQR = abs(eig(A+B*Kd))
%% Testing to see if the controller with gain GP is stable for continuous actions
%Closed_Loop Simulation Using LQR (discrete)
    figure(6)
    [X_hist, U_hist] = simulate(A,B,GP,X0,length(Time));
    subplot(2,1,1);
    plot(Time,Cc*X_hist)
    title(['Closed-Loop MPQL with CONTINUOUS ACTION with r=',num2str(r),' gamma=', num2str(gamma)]);
    
    %Closed-Loop control signal
    subplot(2,1,2);
    plot(Time(1,1:end-1), U_hist)
    title('Control Signal (input) MPQL with CONTINUOUS ACTION');
    
