clear all
%========System 1=========
if (0)
    ac = [0.9385 0.0486; -2.4300 0.9239];
    bc = [0.0012 0.0486]';
    %Discrete
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    %weight matrices
    Q = 2*eye(2);
    R = 1;
    r = 10;          %predictive horizon
    gamma = 0.8;    %discount factor
    x_init = [2;-1];
    inputVals = [-3,3];
end
%=========================
%========System 2=========
if(0)
    m1=1;m2=1;k1=100;k2=100;c1=0.3;c2=-0.3;
    M=[m1 0;0 m2];K=[k1+k2,-k2;-k2,k2];C=[c1+c2,-c2;-c2,c2];Bf=[1 0]';
    ac=[zeros(2,2),eye(2,2);-inv(M)*K,-inv(M)*C];
    bc=[zeros(2,1);inv(M)*Bf];
    c=[1000 0 0 0];d=0;
    dt=0.05;
    [a,b]=c2d(ac,bc,dt);
    Q=[1 0 0 0;
       0 1 0 0;
       0 0 100 0
       0 0 0 100];
    R=1; 
    r = 10;
    gamma=0.8;
end
%=========================
%=====System 3(Marginally Stable)=====
if(1)
    a=[1 0;1 1];b=[1 0]';
    [n,numinputs]=size(b);

    Q=2*eye(n,n);
    R=eye(numinputs,numinputs);
    gamma=1;
end
%======================================
%=====System 4(Marginally stable)======
if(0)
    ac=[0 1; -50 0];
    bc=[0 1]';
    dt = 0.05;
    [a,b] = c2d(ac, bc, dt);
    Q=2*eye(2,2);
    R=1;
    gamma=0.8;
end
%=======================================


if(1) %verifying P
    n = size(a,1);
    numInputs = size(b,2);
    r=2;
%     for j=2:20
%         r = j;
        S = calculateAnalyticalS(a,b,r,gamma,Q,R);

%         Sxu = S(1:n, n+1:n+r*numInputs);
%         Suu = S(n+1:n+r*numInputs,n+1:n+r*numInputs);
%         GS = -pinv(Suu)*Sxu';
%         GL(j,:) = GS(1:numInputs,:);

        %--------P form eq[70]---------
%         P_70 = calculateAnalyticalP(a,b,r,S);
        %----------End of P_70----------

        %---------P form eq[79]---------
        P_79 = calculateNumericalP(a,b,Q,R,r,gamma,S)
        %-------End of P_79-------------
        
        %------- new P matrix -------------
%         P_new = newPmatrix(a,b,r,S);
        %--------end of new P matrix ------
        P = P_79;
        Pxu = P(1:n,n+1:n+numInputs);
        Puu = P(n+1:n+numInputs, n+1:n+numInputs);
        GP = -pinv(Puu)*Pxu'
%         G(j,:) = -pinv(Puu)*Pxu';
%     end

%     plot(diffs(:,2:end))
%     title('Sum of elementwise abs(P70 - P79) vs r');
%     xlabel('r');
%     ylabel('sum(sum(abs(P70 - P79)))');
%     legend('u(k)=GL*x(k)+randn*0.01 and u(k+1)=GL*x(k+1)');
end

% plot(G)
% hold on
% plot(GL)
% legend('G','GL')
% hold off



%--------Sumulation-------
if(0)
    numIter = 200;
    x = zeros(2,numIter);
    u = zeros(1,numIter);
    Qval = zeros(1,size(inputVals, 2));
    x(:,1) = x_init;
    for k=1:numIter
        for i=1:size(inputVals, 2)
           xu = vertcat(x(:,k), inputVals(:,i));
           Qval(:,i) = xu'*P_70*xu;
        end
        [minVal, minIndex] = min(Qval);
        u(:,k) = inputVals(:,minIndex);
        x(:,k+1) = a*x(:,k) + b*u(:,k);
    end
    
    subplot(2,1,1);
    plot(1:numIter+1, x)
    subplot(2,1,2);
    plot(1:numIter, u)
end


