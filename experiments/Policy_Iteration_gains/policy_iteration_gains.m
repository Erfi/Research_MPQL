% This script prodeuces plots for each step of the policy iteration and
% shows how Q learning trough policy itetation can achieve LQR gains.
%--------------------------------------------------------------------------

clear all;
close all;

%--- System ---
[A,B,Cc,Dc,Q,R,Ac,Bc] = getSystemModel(4);
load('trussmodelforErfan.mat','dt');
[n,m] = size(B);
numIters = 6;
r = 40;
gamma = 0.85;
%--------------

S = calculateAnalyticalS(A,B,r,gamma,Q,R);
[~,GL,~] = extractGainFromS(S,n,m);
[P,GP_hist] = calculateOptimalP_PI(A,B,Q,R,r,1.0,GL,false,numIters);
%%  plotting
%---extract first 5 gains---
G_LQR = -dlqr(A,B,Q,R);
G_LQR_extracted = repmat(G_LQR(1,1:5),numIters+1,1);
LQR_plot = plot(0:numIters, G_LQR_extracted,'LineWidth',2.5, 'Color',[0.8,0.5,0.5])
hold on

GP_extracted = GP_hist(1:2:end,1:5);
QL_plot = plot(0:numIters,GP_extracted,'LineWidth',2.5, 'Color',[0.4,0.7,0.4])
hold off

set(gca,'FontSize',25) %set axis properties
xlabel('Iteration Number','FontSize', 30)
ylabel('Gain Values','FontSize', 30)
Leg = legend([LQR_plot(1),QL_plot(1)],'LQR','Q-Learning');
set(Leg, 'FontSize', 20)
grid on;



