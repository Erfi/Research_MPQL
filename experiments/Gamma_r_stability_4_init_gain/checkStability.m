function [ status ] = checkStability( Gain, A, B )
% Helper function to see if a gain is stablizing or not.
%
% Args:
%   Gain: Gain to be examined
%   A: System Dynamic Matrix
%   B: System Input Matrix
%
% Returns:
%   status: 1 if stable, 0 if unstable 
%--------------------------------------------------------------------------
    status = max(abs(eig(A+B*Gain))) < 1;
end

