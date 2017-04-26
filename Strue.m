%Erfan Azad
%

clear all

ac = [0 1; -50 -0.3];
bc = [0 1]';

dt = 0.05;

[a,b] = c2d(ac, bc, dt);
c = eye(2,2);
d = zeros(2,1);

% [n,~] = size(a);
% [numStates, numInputs] = size(b);

r = 3;
Q = 2*eye(2,2);
R = 1;
gamma = 0.5;

S_true = calculateAnalyticalS(a,b,r,gamma,Q,R);

% 
% capQ = zeros(n*r,n*r);
% capGammaQ = zeros(n*r,n*r);
% for i=1:r
%    capQ(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = Q;
%    capGammaQ(n*(i-1)+1:n*i,n*(i-1)+1:n*i)= eye(n,n) * (sqrt(gamma))^(i-1);
% end
% Qgamma = capGammaQ*capQ*capGammaQ;
% 
% capR = zeros(numInputs*r, numInputs*r);
% capGammaR = zeros(numInputs*r, numInputs*r);  
% for i=1:r
%     capR(numInputs*(i-1)+1:numInputs*i, numInputs*(i-1)+1:numInputs*i) = R;
%     capGammaR(numInputs*(i-1)+1:numInputs*i, numInputs*(i-1)+1:numInputs*i) = eye(numInputs, numInputs) * (sqrt(gamma))^(i-1);
% end
% Rgamma = capGammaR*capR*capGammaR;
% 
% % p1 and p2 matrices 
% P1 = zeros(n*r, n);
% for i=1:r
%     P1(n*(i-1)+1:n*i,:) = a^i;
% end
% 
% P2 = zeros(numStates*r, numInputs*r);
% for i=1:r %(for row=i)
%     row = zeros(numStates, numInputs*r);
%     for j=1:i %(for col=j)
%         row(:,numInputs*(j-1)+1:numInputs*j) = a^(i-j)*b;
%     end
%    P2(numStates*(i-1)+1:numStates*i,:) = row;
% end

% 
% %S_gamma
% S_true = [Q+P1'*Qgamma*P1, P1'*Qgamma*P2;
%         P2'*Qgamma*P1, Rgamma+P2'*Qgamma*P2]
