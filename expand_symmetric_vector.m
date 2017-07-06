function [ expanded_matrix ] = expand_symmetric_vector( vector,n)
% Creates a vector of length (l+l^2)/2 and produces a symmetric matrix
% of size (lxl)
%
% Args:
%   vector: Input reduced vector
%   n: Dimention of the output expanded_matrix
%
% Returns:
%   expanded_matrix: A symmetric nxn matrix
%--------------------------------------------------------------------------
expanded_matrix = zeros(n,n);

sparsity_mask = logical(triu(ones(n,n)));
expanded_matrix(sparsity_mask) = vector;

sparsity_mask = logical(triu(ones(n,n),1));
upper_values = expanded_matrix(sparsity_mask);

sparsity_mask = logical(tril(ones(n,n),-1));
expanded_matrix(sparsity_mask) = upper_values;
end

