function [X,U] = simulate(A,B,Gain,X0,numIter)
% This is a funciton that simulates the system given by the A and B matrix
% using the Gain as the controller gain. 
%
% Args:
%   A: System matrix
%   B: Input matrix
%   Gain: Controller gain
%   X0: Initial state of the system
%   numIter: suration of the simulation
%
% Returns:
%   X: Complete state history
%   U: Complete control signal (input) history
%--------------------------------------------------------------------------
    [n,m] = size(B);
    X = zeros(n,numIter);
    U = zeros(m,numIter-1);
    X(:,1) = X0;
    for i=1:numIter-1
       U(:,i) = Gain*X(:,i);
       X(:,i+1) = A*X(:,i) + B*U(:,i);
    end
end

