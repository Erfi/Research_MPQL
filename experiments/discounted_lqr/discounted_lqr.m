% In this experiment we will used the Discrete Alrebraic Raccati Equations
% solver (dare) to compute the discouted-lqr gain at each time step and
% 1) See if it converges
% 2) Does it converge to the same gain as the discounted Q-Learning
%--------------------------------------------------------------------------

clear all;
close all;

%--- System ---
[A,B,Cc,Dc,Q,R,Ac,Bc] = getSystemModel(4);

load('trussmodelforErfan.mat','dt');
[n,m] = size(B);
%--------------

%--- Simulation Variables ---
gamma = 1.0
%---------------------------

%% Discounted LQR
K_t = discounted_dlqr(A,B,Q,R,300,gamma);
%-----Stability------
K_t__stability = checkStability(K_t,A,B)
%% Discounted Q-Learning
r=nan;
G_inital = randn(m,n)*2;
G_initial_stability = checkStability(G_inital,A,B)
P = calculateOptimalP_PI(A,B,Q,R,r,gamma,G_inital,false,10);
G_QL = extractGainFromP(P,n);
%------ Stability -------
G_QL__stability = checkStability(G_QL,A,B)
%% Undiscounted LQR
G_LQR = -dlqr(A,B,Q,R);
%------ Stability --------
G_LQR__stability = checkStability(G_LQR,A,B)
