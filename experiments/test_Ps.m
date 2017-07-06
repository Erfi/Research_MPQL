% This file tests different ways to calculate the P matrix
%--------------------------------------------------------------------------
clear all;
close all;

%--- System Model ---
[a,b,C,D,Q,R,ac,bc] = getSystemModel(3);
%--------------------

[n,m] = size(b);
dt = 0.05;
r = 5;
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

%-------
P_ana_s = calculateAnalyticalPs(a,b,r,S);
GP_ana_s = extractGainFromP(P_ana_s,n)
%-------
P_ana = calculateAnalyticalP(a,b,r,S);
GP_ana = extractGainFromP(P_ana,n)
%-------
P_num_batch = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
GP_num_batch = extractGainFromP(P_num_batch,n)
%-------
P_num_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
GP_num_RLS = extractGainFromP(P_num_RLS,n)
%-------
P_PI = calculateOptimalP_PI(a,b,Q,R,r,gamma,S,true,5);
GP_PI = extractGainFromP(P_PI,n)
%-------
[~,GL,~] = extractGainFromS(S,n,m);
GL
%-------
GLQR = -dlqr(a,b,Q,R)



