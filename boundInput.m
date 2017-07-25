function [ bounded_input ] = boundInput( input, bound )
% Bounds the input such that -bound <= input <= bound
%
% Args:
%   input: Input/Control Signal
%   bound: Bound for the input (Needs to be a positive number)
%
% Returns:
%   bounded_input: bounded input such that -bound <= bounded_input <= bound
%--------------------------------------------------------------------------
    bound = abs(bound); %make sure bound is a positive number
    bounded_input = min(bound, max(input, -bound));
end

