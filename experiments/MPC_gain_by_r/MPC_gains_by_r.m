% This script is for plotting MPC gains for different predictive horizons
% and compare them to the LQR gains.
%--------------------------------------------------------------------------

clear all;
close all;

%--- System ---
[A,B,Cc,Dc,Q,R,Ac,Bc] = getSystemModel(4);

load('trussmodelforErfan.mat','dt');
[n,m] = size(B);
%--------------

%--- Simulation Variables ---
gamma = 1.0;
r_vals = 10:10:250;
% only the forst five gains
MPC_gains = zeros(5,length(r_vals));
%---------------------------

%------- MPC gains --------
index = 1;
for r=r_vals
    G_MPC = getMPCGain(A,B,Q,R,r);
    MPC_gains(:,index) = G_MPC(1,1:5);
    index = index+1;
end
%---------------------------

%------- LQR gains ---------
G_LQR = -dlqr(A,B,Q,R);
LQR_gains = G_LQR(1,1:5);
%---------------------------

%% plots

LQR_plot = plot(r_vals, repmat(LQR_gains', 1,length(r_vals)), 'LineWidth',2.5, 'Color',[0.8,0.5,0.5]);
hold on;
MPC_plot = plot(r_vals, MPC_gains','LineWidth',2.5, 'Color',[0.3,0.3,0.5]);
hold off;
set(gca,'FontSize',25) %set axis properties
xlabel('Prediction Horizon','FontSize', 30)
ylabel('Gain Values','FontSize', 30)
Leg = legend([LQR_plot(1),MPC_plot(1)],'LQR','MPC');
set(Leg, 'FontSize', 20)
grid on;