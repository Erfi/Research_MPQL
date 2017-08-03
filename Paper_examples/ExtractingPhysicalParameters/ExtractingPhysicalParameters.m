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
runDiscreteLQR =   false;
runMPC_Original =  false;
runMPC =           false;
runImplicitMPQL =  true;
runImplicitMPQL_Continuous = false;
%---------------------------------------------------

%------------Simulation Variables-------------------
% Time = 0:dt:2;
Time = 0:500;
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
title('Open-Loop Simulation');
%Closed_Loop Simulation Using LQR (discrete)
Kd = -dlqr(A,B,Q,R)
[X_hist_LQR_D, U_hist_LQR_D] = simulate(A,B,Kd,X0,length(Time));
subplot(3,1,2);
plot(Time,Cc*X_hist_LQR_D);
title('Closed-Loop Simulation LQR (discrete-time)');

%Closed-Loop control signal
subplot(3,1,3);
plot(Time(1,1:end-1), U_hist_LQR_D);
title('Control Signal (input) LQR (discrete-time)');
end
%------------------------------------------------------

%------------------MPC (Original)----------------------
if (runMPC_Original)
    r = 40;
    G_MPC = getMPCGain(A,B,Q,R,r)

    %Open_Loop Simulation Using LQR (discrete)
    [Y,~]=dlsim(A,B,Cc,Dc,U,X0);
    figure(3);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPC\_Original');

    %Closed_Loop Simulation Using LQR (discrete)
    [X_hist_MPC_Orig, U_hist_MPC_Orig] = simulate(A,B,G_MPC,X0,length(Time));
    subplot(3,1,2);
    plot(Time,Cc*X_hist_MPC_Orig)
    title(['Closed-Loop simulation MPC\_Original with r=',num2str(r)]);
    
    %Closed-Loop control signal
    subplot(3,1,3);
    plot(Time(1,1:end-1), U_hist_MPC_Orig)
    title('Control Signal (input) MPC\_Original');
end

%------------------------------------------------------

%------------------MPC (from S)------------------------
if(runMPC)
    r = 40;
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
    [X_hist_MPC_S, U_hist_MPC_S] = simulate(A,B,GL,X0,length(Time));
    subplot(3,1,2);
    plot(Time,Cc*X_hist_MPC_S)
    title(['Closed-Loop simulation MPC (discrete-time) with r=',num2str(r),' gamma=', num2str(gamma)]);
    
    %Closed-Loop control signal
    subplot(3,1,3);
    plot(Time(1,1:end-1), U_hist_MPC_S)
    title('Control Signal (input) MPC (discrete-time)');
end
%--------------------------------------------------------

%--------------------MPQL(implicit)----------------------
if(runImplicitMPQL)
    
    r = 40;
    gamma = 0.85;
    
  
    input_vals = [-500000:100000:500000,-10000:1000:10000,-1000:100:1000,-100:10:100];
    input_vals = [-40000,40000,-30000,30000,-20000,20000,-10000,10000,-1000,1000,-100,100,0];
    input_vals = [-5000, 5000, -500,500, -50, 50, 0];
    input_vals = [-5,5,-3,3,-1,1,-0.5,0.5,0];
    input_vals = repmat(input_vals, [m,1]);
        
    numIter = length(Time);
    [X_hist_MPQL, U_hist_MPQL, GP,P] = implicitMPQL(A,B,Q,R,r,gamma,X0,input_vals, numIter,1, false);
    GP
    
%     Open_Loop Simulation
    [Y,X]=dlsim(A,B,Cc,Dc,U,X0);
    figure(5);
    subplot(3,1,1);
    plot(Time,Y);
    title('Open-Loop simulation MPQL (discrete-time)'); 
    
    
%     Closed_Loop Simulation
    subplot(3,1,2);
    plot(Time,Cc*X_hist_MPQL(:,1:end-1))
    title(['Closed-Loop simulation MPQL (discrete-time) with r=', num2str(r),' gamma=', num2str(gamma)]);
    %Control signal
    subplot(3,1,3);
    plot(Time, U_hist_MPQL);
    title('Control Input');


%     subplot(2,1,1);
%     plot(Time,Cc*X_hist_MPQL(:,1:end-1))
%     set(gca,'FontSize',25) %set axis properties
%     xlabel('Time Step','FontSize', 30)
%     ylabel('State', 'FontSize', 30);
%     %Control signal
%     subplot(2,1,2);
%     plot(Time, U_hist_MPQL);
%     set(gca,'FontSize',25) %set axis properties
%     xlabel('Time Step','FontSize', 30)
%     ylabel('Control Input', 'FontSize', 30);

end
%---------------------------------------------------------

%--------MPQL(implicit) with coninuous input--------------
if(runImplicitMPQL_Continuous)
    figure(6)
    [X_hist_MPQL_C, U_hist_MPQL_C] = simulate(A,B,GP,X0,length(Time));
    subplot(2,1,1);
    plot(Time,Cc*X_hist_MPQL_C)
    title(['Closed-Loop MPQL with CONTINUOUS ACTION with r=',num2str(r),' gamma=', num2str(gamma)]);
    
%     set(gca,'FontSize',25) %set axis properties
%     xlabel('Time Step','FontSize', 30)
%     ylabel('State', 'FontSize', 30); 
    
    %Closed-Loop control signal
    subplot(2,1,2);
    plot(Time(1,1:end-1), U_hist_MPQL_C)
    title('Control Signal (input) MPQL with CONTINUOUS ACTION'); 
    
%     set(gca,'FontSize',25) %set axis properties
%     xlabel('Time Step','FontSize', 30)
%     ylabel('Control Input', 'FontSize', 30);
end
%----------------------------------------------------------

%% Looking at the magnitude of the eigenValues of GP:
format short
eigMag = abs(eig(A+B*GP))
% eigMagLQR = abs(eig(A+B*Kd))
%% -----------Creating plots for paper--------------
subplot(3,2,[1,2])
%Open-loop simulation
[Y,X]=dlsim(A,B,Cc,Dc,U,X0);
plot(Time,Y);
title('Open-Loop simulation'); 
xlabel('Time (s)');
ylabel('State');

%LQR
%Closed_Loop
subplot(3,2,3);
plot(Time,Cc*X_hist_LQR_D);
title('LQR');
ylabel('State');
%control signal
subplot(3,2,5);
plot(Time(1,1:end-1), U_hist_LQR_D);
ylabel('Input');
xlabel('Time (s)')

%Q-Learning
%Closed_Loop
subplot(3,2,4)
plot(Time,Cc*X_hist_MPQL_C)
title(['Q-Learning with gamma=',num2str(gamma)]);
ylabel('State')
%control signal
subplot(3,2,6);
plot(Time(1,1:end-1), U_hist_MPQL_C)
ylabel('Input')
xlabel('Time (s)')

%% Visualizing CTG and Utility

if (0)
%--- LQR ----
[~, Utility_LQR] = getCostToGo(X_hist_LQR_D,U_hist_LQR_D,false ,Q , R);
figure;
plot(Time(1:size(Utility_LQR,2)), Utility_LQR, 'LineWidth',3)
hold on
area(Time(1:size(Utility_LQR,2)), Utility_LQR,'FaceAlpha', 1 , 'EdgeColor', [0.2,0.2,0.6] ,'FaceColor', [0.7,0.15,0.15])
hold off
grid on;
title('U(k) over time for LQR')
xlabel('Time(s)')
ylabel('U(k)')

Area_LQR = sum(Utility_LQR);
disp(['V(k) or Area under the curve for LQR: ',num2str(Area_LQR)]);
%-------------
end

if (0)
%--- MPC Original ----
[~, Utility_MPC] = getCostToGo(X_hist_MPC_Orig, U_hist_MPC_Orig, false, Q, R);
r = 50;
figure;
plot(Time(1:size(Utility_MPC,2)), Utility_MPC, 'LineWidth',3)
hold on
area(Time(1:size(Utility_MPC,2)), Utility_MPC,'FaceAlpha', 1 , 'EdgeColor', [0.2,0.2,0.6] ,'FaceColor', [0.7,0.15,0.15])
hold off
axis([0, 2, 0, 2*1e4])
grid on;
title(['U(k) over time for MPC with r=',num2str(r)])
xlabel('Time(s)')
ylabel('U(k)')

Area_MPC = sum(Utility_MPC);
disp(['V(k) or Area under the curve for MPC: ',num2str(Area_MPC)]);
%-------------
end


if (1)
%--- MPC from S ----
[~, Utility_MPC_S] = getCostToGo(X_hist_MPC_S, U_hist_MPC_S, false, Q, R);
figure;
plot(Time(1:size(Utility_MPC_S,2)), Utility_MPC_S, 'LineWidth',3)
hold on
area(Time(1:size(Utility_MPC_S,2)), Utility_MPC_S,'FaceAlpha', 1 , 'EdgeColor', [0.2,0.2,0.6] ,'FaceColor', [0.7,0.15,0.15])
hold off
axis([0, 500, 0, 2*1e4])

Area_MPC_S = sum(Utility_MPC_S);
disp(['V(k) or Area under the curve for MPC_S: ',num2str(Area_MPC_S)]);

grid on;
% title('Utility for MPC with r=40 and gamma=1.0')
set(gca,'FontSize',25) %set axis properties
xlabel('Time Step','FontSize', 30)
ylabel('U(k)','FontSize', 30)

%-------------
end


if (0)
%--- MPQL ---
[~, Utility_MPQL] = getCostToGo(X_hist_MPQL, U_hist_MPQL, false, Q, R);
figure;
plot(Time(1:size(Utility_MPQL,2)), Utility_MPQL, 'LineWidth',3)
hold on
area(Time(1:size(Utility_MPQL,2)), Utility_MPQL,'FaceAlpha', 1 , 'EdgeColor', [0.2,0.2,0.6] ,'FaceColor', [0.7,0.15,0.15])
hold off
axis([0, 500, 0, 2*1e4])
grid on;
title(['U(k) over time for MPQL',' iteration 5'])
set(gca,'FontSize',25) %set axis properties
xlabel('Time Step','FontSize', 30)
ylabel('U(k)','FontSize', 30)

Area_MPQL = sum(Utility_MPQL);
disp(['V(k) or Area under the curve for MPQL: ',num2str(Area_MPQL)]);
%-------------
end