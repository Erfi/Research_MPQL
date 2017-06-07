% This file tests different ways to calculate the S matrix
%--------------------------------------------------------------------------
clear all;
close all;

% 2 DOF system
if(1)
    disp('2 DOF system is selected')
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
end

% 8 DOF system
if(0)
    disp(' 8 DOF system is selected')
    n = 8; %number of states (degrees of freedom)
    m = 1; %number of inputs (we are assuming single input)
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
    Bf(1,1) = 1; %assuming single input
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








