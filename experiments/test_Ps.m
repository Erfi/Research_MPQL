% This file tests different ways to calculate the P matrix
%--------------------------------------------------------------------------
clear all;
close all;

%--- System Model ---
% [a,b,C,D,Q,R,ac,bc] = getSystemModel(3);
%--------------------
[a,b,C,D,Q,R,ac,ac] = getSystemModel(4);
load('trussmodelforErfan.mat','dt');

[n,m] = size(b);
dt = 0.004;
r = 10;
gamma = 1;

%--- Checking for controlibility ---
Co = ctrb(a,b);
if rank(Co) == n 
    disp('System is controlable')
else
    disp('System is not controllable')
end

format short
[U,S,V] = svd(Co);

[maxVal, Imax] = max(diag(S));
singleVals = diag(S)/maxVal;
semilogy(singleVals,'*')
title('Normalized singular values of the Controlability matrix')
xlabel('Number of singular values')
ylabel('Value of singular values')
%-----------------------------------


S = calculateNumericalS(a,b,r,gamma,Q,R);
%%
%-------
[~,GL,~] = extractGainFromS(S,n,m);
GL
%-------
P_ana_s = calculateAnalyticalPs(a,b,r,S);
GP_ana_s = extractGainFromP(P_ana_s,n)
%-------
P_ana = calculateAnalyticalP(a,b,r,S);
GP_ana = extractGainFromP(P_ana,n)
%-------
P_num_batch = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
GP_num_batch = extractGainFromP(P_num_batch,n)
% abs(eig(a+b*GP_num_batch))
%-------
P_num_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
GP_num_RLS = extractGainFromP(P_num_RLS,n)
%-------
% My_Arbitrary_Small_Gain = randn(size(GL))*1e-3;
[P_PI,GPs_PI] = calculateOptimalP_PI(a,b,Q,R,r,gamma,S,true,5);
GP_PI = extractGainFromP(P_PI,n)
% abs(eig(a+b*GP_PI))
%-------
GLQR = -dlqr(a,b,Q,R)
% abs(eig(a+b*GLQR))



