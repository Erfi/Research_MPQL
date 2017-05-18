function [X1] = simStep(A,B,U0,X0)
% This is a funciton that simulates the system for one step given A and B 
% using the Gain as the controller gain. 
%
% Args:
%   A: System matrix
%   B: Input matrix
%   Gain: Controller gain
%   X0: Initial state of the system
%   numIter: duration of the simulation
%
% Returns:
%   X1: Next State
%--------------------------------------------------------------------------
    X1 = A*X0 + B*U0;
end

