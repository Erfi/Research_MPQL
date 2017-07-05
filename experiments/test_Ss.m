% This file tests different ways to calculate the S matrix
%--------------------------------------------------------------------------
clear all;
close all;

[a,b,C,D,Q,R,ac,bc] = getSystemModel(3);

[n,m] = size(b);
r = 10;
gamma = 1;
dt = 0.05;

%-------
S_ana = calculateAnalyticalS(a,b,r,gamma,Q,R);
Sxu = S_ana(1:n, n+1:n+r*m);
Suu = S_ana(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL_ana = G(1:m,:) 
%-------
S_num_batch = calculateNumericalS(a,b,r,gamma,Q,R);
Sxu = S_num_batch(1:n, n+1:n+r*m);
Suu = S_num_batch(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL_num_batch = G(1:m,:) 
%--------
S_num_RLS = calculateNumericalS_RLS(a,b,r,gamma,Q,R);
Sxu = S_num_RLS(1:n, n+1:n+r*m);
Suu = S_num_RLS(n+1:n+r*m,n+1:n+r*m);
G = -pinv(Suu)*Sxu';
GL_num_RLS = G(1:m,:) 
%--------
GLQR = -dlqr(a,b,Q,R)







