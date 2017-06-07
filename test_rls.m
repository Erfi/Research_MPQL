% script for testing the RLS algorithm
%--------------------------------------------------------------------------
m = 1000; %number of rows (data points)
n = 2;    %number of cols (features)
X = rand(m,n); % same as LHS
theta = rand(n,1) %estimates (P matrix in vertical form)
Y = X*theta; % data that fits the A and theta

theta_estimate = zeros(n,m+1); %to hold all the estimates
P = eye(n,n); %covarience matrix inv(X'X)

for k=1:m
   [P, theta_estimate(:,k+1)] = rls_one_step(P,theta_estimate(:,k),X(k,:), Y(k,:)); 
end

plot(theta_estimate')