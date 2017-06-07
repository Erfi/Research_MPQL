function [ P, theta ] = rls_one_step( P,theta,x,y )
% Calculates one step of the RLS algorithm using covarience matrix.
%
% Args:
%   P:     Covarience at time k
%   theta: The estimate at tiem k
%   x:     New system row
%   y:     New data row
%
% Returns:
%   P:     Covarience at time k+1
%   theta: Estimate at time k+1
%--------------------------------------------------------------------------
K = P*x'./(1+x*P*x'); %gain
P = P-K*(x*P); %update covarience matrix
theta = theta + K*(y-x*theta); %update estimate
end

