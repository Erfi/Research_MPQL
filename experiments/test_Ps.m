% This file tests different ways to calculate the P matrix
%--------------------------------------------------------------------------
clear all;
close all;

%--- System Model ---
[a,b,C,D,Q,R,ac,bc] = getSystemModel(2);
%--------------------

[n,m] = size(b);
dt = 0.05;
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
%-----------------------------------


S = calculateAnalyticalS(a,b,r,gamma,Q,R);

%-------
P_ana_s = calculateAnalyticalPs(a,b,r,S);
Pxu = P_ana_s(1:n, n+1:end);
Puu = P_ana_s(n+1:end, n+1:end);
GP_ana_s = -pinv(Puu)*Pxu'
%-------
P_ana = calculateAnalyticalP(a,b,r,S);
Pxu = P_ana(1:n, n+1:end);
Puu = P_ana(n+1:end, n+1:end);
GP_ana = -pinv(Puu)*Pxu'
%-------
P_num_batch = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
Pxu = P_num_batch(1:n, n+1:end);
Puu = P_num_batch(n+1:end, n+1:end);
GP_num_batch = -pinv(Puu)*Pxu'
%-------
P_num_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
Pxu = P_num_RLS(1:n, n+1:end);
Puu = P_num_RLS(n+1:end, n+1:end);
GP_num_RLS = -pinv(Puu)*Pxu'
%-------
P_PI = calculateOptimalP_PI(a,b,Q,R,r,gamma,S,true,5);
Pxu = P_num_RLS(1:n, n+1:end);
Puu = P_num_RLS(n+1:end, n+1:end);
GP_PI = -pinv(Puu)*Pxu'
%-------
G = extractGainFromS(S,n,m);
GL = G(1:m,:)
%-------
GLQR = -dlqr(a,b,Q,R)



