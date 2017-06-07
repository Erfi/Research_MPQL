function [ A,B,C,D,Q,R ] = getSystemModel( sysNumber )
% Since writing out the system model in every driver file can makes things
% messy this function will give you the system model that you request. 
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
%---------------------------------
%----- System #1 ----
% A 2 DOF system
%--------------------
%
%----- System #2 ----
% An 8 DOF system, 4 inputs
%--------------------
%--------------------------------------------------------------------------

if (sysNumber == 1)
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [A,B] = c2d(ac, bc, dt);
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
                  -17.7 -17.7 0 0 0 -50.0 0 86.4]      * 0.0; % decrease damping

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
    C = eye(1,2*n); % n-output 
    D = 0; % direct transition matrix
    Q = eye(2*n); 
    R = eye(m)*1e-4;
    dt = 0.05; % sampling delta using ~6 * highest frequency of the system
    [A,B] = c2d(Ac, Bc, dt); % discrete system dynamic
else
    disp('Not a system number. type "help getSystemModel" for more info')
end
end

