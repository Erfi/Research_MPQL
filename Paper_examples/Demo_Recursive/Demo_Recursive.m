% This is a demonstration of recursively finding S --> GL --> P --> repeat 
% to control a system in real time.
clear all;
%------System (marginally stable)------
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 5;
    gamma=1.0; %1.0 --> LQR uses gamma = 1.0
%------------------

%----Calculating S using different methods-----
% S_ana = calculateAnalyticalS(a,b,r,gamma,Q,R);
% S_num_batch = calculateNumericalS(a,b,r,gamma,Q,R);
% S_num_RLS = calculateNumericalS_RLS(a,b,r,gamma,Q,R);
%----Calculating P using different methods-----
% P_batch = calculateNumericalP(a,b,Q,R,r,gamma,S_ana);
% P_RLS = calculateNumericalP_RLS(a,b,Q,R,r,gamma,S_num_RLS, true);

%--Demo--
[n,m] = size(b);
S = calculateNumericalS_RLS(a,b,r,gamma,Q,R); %using small r
% GL from S
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
GS = -pinv(Suu)*Sxu';
GL = GS(1:m,:);
GP = GL
numIter =10;
for i = 1:numIter
    P = calculateNumericalP_RLS(a,b,Q,R,r,gamma,GP(i,:),false);
    P_xu = P(1:n, n+1:end);
    P_uu = P(n+1:end, n+1:end);
    GP(i+1,:) = -inv(P_uu)*P_xu'
end

%LQR gain
GLQR = -dlqr(a,b,Q,R);

%plotting
plot(1:numIter+1,GP)
hold on;
plot(1:numIter+1,[ones(numIter+1,1)*GLQR(1), ones(numIter+1,1)*GLQR(2)])
hold off;
xlabel('Iteration Number');
ylabel('Gain Values');
title(['LQR gain vs. Recursive GP gain using r=', num2str(r)]);
legend('RLS_1','RLS_2','LQR_1','LQR_2');

