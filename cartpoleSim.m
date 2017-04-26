% Erfan Azad <erfan@dartmouth.edu>
% Date: 3 April 2017
%---------------------------------
clear all;
close all;

m = 0.1; %kg
M = 1;   %kg
l = 1;   %m
g = 9.8; %m/s^2


ac = [0 1 0 0;
    0 0 -m*g/M 0;
    0 0 0 1;
    0 0 g*(m+M)/(M*l) 0];
bc = [0 1/M 0 -1/(M*l)]';
dt = 0.05;
[a,b] = c2d(ac,bc,dt);
Q = eye(4);
R = 1;
r = 30;
gamma = 1;

% counter = 1;
% for r = 20:1:25
%     for gamma = 0.6:0.05:1
%         r
%         gamma
        S = calculateNumericalS(a,b,r, gamma, Q,R);
        P = calculateNumericalP(a,b,Q,R,r,gamma,S);
%         Pnew = newPmatrix(a,b,r,S);
        pos = 1;
        pos_dot = 0; 
        theta = 1;
        theta_dot = 0;
        state = [pos, pos_dot, theta, theta_dot]';
        u_options = [-15:1:15];
        
        %-----GL------
        n = size(a,1);
        numInputs = size(b,2);
        Sxu = S(1:n, n+1:n+r*numInputs);
        Suu = S(n+1:n+r*numInputs,n+1:n+r*numInputs);
        G = -pinv(Suu)*Sxu';
        GL = G(1:numInputs,:);
        %-------------

        for i=1:400
            for j=1:size(u_options,2)
                xu = vertcat(state(:,i),u_options(j));
                Qval(j) = xu'*P*xu;
            end
            [minVal, minIndex] = min(Qval);
            u_optimal(i) = u_options(minIndex);
%             if(i>224)
%                 u_optimal(i) = 0;
%             end
%             u_optimal(i) = GL*state(:,i); %Continuous
            state(:,i+1) = a*state(:,i) + b*u_optimal(i);
        end
%         if(sum(abs(state(:,end))) > 2)
%             result(:,counter) = [r gamma -1]';
%         else
%             result(:,counter) = [r gamma 1]';
%         end
%         counter = counter + 1;
%     end
% end

% scatter(result(1,:), result(2,:), 50, result(3,:), 'filled')
% title('Stabilization success or failure for different values of r and gamma')
% xlabel('r')
% ylabel('gamma')
% legend('Failure')

CTG = 0;
for k = 1:400
   CTG = CTG + state(:,i+1)'*Q*state(:,i+1) + u_optimal(i)'*R*u_optimal(i); 
end

if(1)
%     state(3:3,:) = rad2deg(state(3,:));
    subplot(2,1,1)
    plot(state')
    title(['cart-pole stabilization r=',num2str(r),' gamma=',num2str(gamma)])
    xlabel('time step')
    legend('z','z-dot', 'theta', 'theta-dot')
    subplot(2,1,2)
    plot(u_optimal)
    title(['control input. u-options[' , num2str(u_options), ']'])
    xlabel('time step')
    ylabel('control value')
    legend('input')
end