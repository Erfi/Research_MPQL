function [ f_combined ] = combineSymmetricCoefficientsInRow( f )
% This function takes in an array f with size (n^2) which
% will be later multiplied with a symetric (n*n) matrix
% and adds the coefficients corresponding to the symmetric
% elements of the (n*n) matrix.
% 
% Args:
%   f: input array of size (n^2)
% 
% Returns:
%   f_combined: an array of size (((n^2)-n)/2)+n
%
    n = sqrt(size(f,2));
    f_matrix = reshape(f, [n,n]);
    f_combined = zeros(1,((size(f,2)-n)/2)+n);
    index = 1;
    for i=1:n
        for j=i:n
            if (i==j)
                f_combined(1,index) = f_matrix(i,j);
            else
                f_combined(1,index) = f_matrix(i,j) + f_matrix(j,i);
            end
            index = index + 1;
        end
    end
end

