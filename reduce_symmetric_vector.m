function [ reduced_vector, sparsity_mask ] = reduce_symmetric_vector( vector )
% Takes in an array that will be multiplied to a symmetric matrix P. Only
% keeps the elements corresponding to the upper triangle of P. And
% doubles the value of all the symmetric elements.
%
% Args:
%   vector: An array/vector containing n elements (1xn)
%   
% Returns:
%   reduced_vector: An array/Vector containing k=(n+n^2)/2 elements (1xk)
%--------------------------------------------------------------------------

len = max(size(vector));
dim = sqrt(len);
reduced_len = (dim+(dim^2))/2;

matrix = reshape(vector,[dim,dim]);
upper = triu(matrix,1)*2; %double the value of the elements above diagnal
diagonal = eye(dim).*matrix;
new_matrix = upper+diagonal;

sparsity_mask = logical(triu(new_matrix));
reduced_vector = new_matrix(sparsity_mask);
reduced_vector = reshape(reduced_vector, [1,reduced_len]);
end

