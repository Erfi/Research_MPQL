% This file tests different ways to calculate the P matrix
%--------------------------------------------------------------------------
clear all;
close all;

% 2 DOF system
if(0)
    disp('2 DOF system is selected')
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
end

% 8 DOF system
if(1)
    disp(' 8 DOF system is selected')
    n = 8; %number of states (degrees of freedom)
    m = 4; %number of inputs (we are assuming single input)
    massMatrix = eye(8)*100;
    stiffnessMatrix = [27071.1 0 0 0 -10000.0 0 -3535.5 -3535.5;
                       0 17071.1 0 -10000.0 0 0 -3535.5 -3535.5;
                       0 0 27071.1 0 -3535.5 3535.5 -10000.0 0;
                       0 -10000.0 0 17071.1 3535.5 -3535.5 0 0;
                       -10000.0 0 -3535.5 3535.5 27071.1 0 0 0;
                       0 0 3535.5 -3535.5 0 17071.1 0 -10000.0;
                       -3535.5 -3535.5 -10000.0 0 0 0 27071.1 0;
                       -3535.5 -3535.5 0 0 0 -10000.0 0 17071.1];
     dampingMatrix = [136.4 0 0 0 -50.0 0 -17.7 -17.7;
                      0 86.4 0 -50.0 0 0 -17.7 -17.7;
                      0 0 136.4 0 -17.7 17.7 -50.0 0;
                      0 -50.0 0 86.4 17.7 -17.7 0 0;
                      -50.0 0 -17.7 17.7 136.4 0 0 0;
                      0 0 17.7 -17.7 0 86.4 0 -50.0;
                      -17.7 -17.7 -50.0 0 0 0 136.4 0;
                      -17.7 -17.7 0 0 0 -50.0 0 86.4];
     dampingMatrix = zeros(n,n);

    %---System Dynamic---
    Ac = [zeros(n,n), eye(n);
          -inv(massMatrix)*stiffnessMatrix -inv(massMatrix)*dampingMatrix];
    Bf =  zeros(n,m);
    Bf(1,1) = 1;
    Bf(2,2) = 1;
    Bf(7,3) = 1;
    Bf(8,4) = 1;
    
    Bc = [zeros(n,m);
          inv(massMatrix)*Bf];
    Cc = eye(1,2*n); % n-output 
    Dc = 0; % direct transition matrix
    Q = eye(2*n); 
    R = 1*1e-4;
    dt = 0.05; % sampling delta using ~6 * highest frequency of the system
    [a,b] = c2d(Ac, Bc, dt); % discrete system dynamic
    end
%==========================================================================
[n,m] = size(b);
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
GP_ana_s = -inv(Puu)*Pxu'
%-------
P_ana = calculateAnalyticalP(a,b,r,S);
Pxu = P_ana(1:n, n+1:end);
Puu = P_ana(n+1:end, n+1:end);
GP_ana = -inv(Puu)*Pxu'
%-------
P_num_batch = calculateNumericalP(a,b,Q,R,r,gamma,S,true);
Pxu = P_num_batch(1:n, n+1:end);
Puu = P_num_batch(n+1:end, n+1:end);
GP_num_batch = -inv(Puu)*Pxu'
%-------
P_num_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S,true);
Pxu = P_num_RLS(1:n, n+1:end);
Puu = P_num_RLS(n+1:end, n+1:end);
GP_num_RLS = -inv(Puu)*Pxu'


