% we have noticed that the Value Iteration scheme does not
% converge. Here we are trying to figure out why. And for that purpose we
% are checking the rank of the LHS matrix.
%------------------------------------------------------------------------
clear all;
close all;

%------System (marginally stable)------
if(1)
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=eye(2);
    R=1;
    r = 5;
    gamma=1; %1.0 --> LQR uses gamma = 1.0
end
%------------------
[n,m] = size(b);
%------------------
S = calculateNumericalS_RLS(a,b,r,gamma,Q,R); %using small r
% GL from S
Sxu = S(1:n, n+1:n+r*m);
Suu = S(n+1:n+r*m,n+1:n+r*m);
GS = -pinv(Suu)*Sxu';
GL = GS(1:m,:) % Initial stable policy/gain
%------------------
G_LQR = -dlqr(a,b,Q,R)
%------------------
[P_PI,GP_PI] = calculateOptimalP_PI(a,b,Q,R,r,gamma,GL,10);
[P_VI,GP_VI, X, U, LHS, RHS] = calculateOptimalP_VI(a,b,Q,R,r,gamma,GL);
GP_VI_final = GP_VI(end,:)

%----closed loop eig values of all gains----
for i = 1:length(GP_VI)
    try
        eigs(i,:) = eig(a+b*GP_VI(i,:));
    catch
        warning('CP_VI has NaN values in the end');
    end
end
%% ----plots-----
stepsBeforeEnding = 5;

figure(1);
subplot(2,1,1);
plot(X(1:end-stepsBeforeEnding,:));
title('x history ::: Terminated few steps before ending');
xlabel('time steps');
ylabel('X');
subplot(2,1,2);
plot(U(1:end-stepsBeforeEnding,:));
title('u history ::: Terminated few steps before ending');
xlabel('time steps');
ylabel('Control Signal');



%plotting
numIter = length(GP_VI)-stepsBeforeEnding;
figure;
plot(1:numIter,GP_VI(1:end-stepsBeforeEnding,:))
hold on;
plot(1:numIter+1,[ones(numIter+1,1)*G_LQR(1), ones(numIter+1,1)*G_LQR(2)])
hold off;
xlabel('Iteration Number');
ylabel('Gain Values');
title(['LQR gain vs. Recursive GP gain using r=', num2str(r)]);
legend('GP-VI_1','GP-VI_2','LQR_1','LQR_2');

%% ---Discussion---
% We can see that when the unstable gains are skipped, the rank of the LHS 
% matrix becomes full. And the gain value converge although not necessarily
% to the LQR gain.