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
gamma = 0.9
%---------------------------

%% Discounted LQR
N = 300;
P_tp1 = (gamma^(N+1))*Q;
for t=N:-1:0
    K_t = -pinv((gamma^(t))*R + B'*P_tp1*B)*B'*P_tp1*A;
    P_t = (gamma^(t))*Q + A'*P_tp1*A - A'*P_tp1*B*pinv((gamma^(t))*R + B'*P_tp1*B)*B'*P_tp1*A;
    P_tp1 = P_t;
end
%-----Stability------
K_t__stability = checkStability(K_t,A,B)
%% Discounted Q-Learning
r=nan;
G_inital = randn(m,n)*1e-3;
G_initial_stability = checkStability(G_inital,A,B)
P = calculateOptimalP_PI(A,B,Q,R,r,gamma,G_inital,false,10);
G_QL = extractGainFromP(P,n);
%------ Stability -------
G_QL__stability = checkStability(G_QL,A,B)
%% Undiscounted LQR
G_LQR = -dlqr(A,B,Q,R);
%------ Stability --------
G_LQR__stability = checkStability(G_LQR,A,B)
