function [ A,B,C,D,Q,R,Ac,Bc ] = getSystemModel( sysNumber )
% Since writing out the system model in every driver file can make things
% messy, this function will give you the system model that you request. 
%
% Args:
%   sysNumber: Number of your desired system
%
% Returns:
%   A: System Dynamic Matrix (Discrete)
%   B: System Input Matrix (Discrete)
%   C: Output Matrix
%   D: Direct Transition Matrix
%   Q: System Weight Matrix
%   R: input Weight Matrix
%   Ac: System Dynamic Matrix (Continuous)
%   Bc: System Input Matrix (Continuous)
%---------------------------------
%----- System #1 ----
% A 1 DOF system, 1 input
%--------------------
%
%----- System #2 ----
% An 8 DOF system, 4 input
%--------------------
%
%----- System #3 ----
% A 2 DOF system, 2 input
%--------------------
%--------------------------------------------------------------------------

if (sysNumber == 1)
    Ac=[0 1; -50 0];
    Bc=[0 1]';
    dt = 0.05;
    [A,B] = c2d(Ac, Bc, dt);
    C = nan;
    D = nan;
    Q=eye(2);
    R=1;
    
elseif (sysNumber == 2)
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
                  -17.7 -17.7 0 0 0 -50.0 0 86.4] * 0.0; %decreasing the damping

    %---System Dynamic---
    Ac = [zeros(n,n), eye(n);
      -inv(massMatrix)*stiffnessMatrix -inv(massMatrix)*dampingMatrix];
    Bf =  zeros(n,m);
    Bf(1,1) = 1;
    Bf(2,2) = 1;
    Bf(7,3) = 1;
    Bf(8,4) = 1;
    Bc = [zeros(n,m);
      pinv(massMatrix)*Bf];
    C=eye(n,2*n);
    D = zeros(size(C,1),m); % direct transition matrix
    Q = eye(2*n); 
    R = 1e-5*eye(m);
    dt = 0.05; % sampling delta using ~6 * highest frequency of the system
    [A,B] = c2d(Ac, Bc, dt); % discrete system dynamic

elseif (sysNumber == 3)
    n = 2; % DOF
    m = 2; % number of inputs
    
    m1 = 1;
    m2 = 1;
    k1 = 100;
    k2 = 100;
    c1 = 0.1;
    c2 = 0.1;
    Mass = [m1 0;
            0 m2];
    Striffness = [k1+k2 -k2;
                  -k2 k2];
    Damping = [c1+c2 -c2;
                -c2 c2];
    
    Ac = [zeros(n,n), eye(n);
      -inv(Mass)*Striffness -inv(Mass)*Damping];
    Bf =  zeros(n,m);
    Bf(1,1) = 1;
    Bf(2,2) = 1;
    Bc = [zeros(n,m);
      pinv(Mass)*Bf];
    
    C = eye(2*n); % all-output 
    D = zeros(size(C,1), m); % direct transition matrix
    dt = 0.05;
    [A,B] = c2d(Ac, Bc, dt); % discrete system dynamic
    Q = eye(2*n)*10;
    R = eye(m);
else
    disp('Not a system number. type "help getSystemModel" for more info')
end
end

