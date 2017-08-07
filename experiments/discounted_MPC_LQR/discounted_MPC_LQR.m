% This is to compare the gains from discounted LQR and discounted MPC.
% We are expecting the gains to be the same for large MPC r.
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

%--- Discounted LQR ---
G_LQR = discounted_dlqr(A,B,Q,R,300,gamma)
%----------------------

%--- Discounted MPC ---
G_MPC = getMPCGain(A,B,Q,R,40,gamma)
%----------------------

%--- Discounted MPC_S ---
% S = calculateAnalyticalS(A,B,250,gamma,Q,R);
% [~,G_MPC_S,~] = extractGainFromS(S,n,m)
%------------------------