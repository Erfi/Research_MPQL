function [ m_combined ] = combineSymmetricCoefficientsInMatrix( m )
% This function takes in a k*(n^2) size matrix m and
% combines symmetric coefficients in each row corresponding
% to an (n*n) symmetric matrix.
%
% Args:
%   m: k*(n^2) matrix of coefficients
%
% Returns:
%   m_combined: k*((((n^2)-n)/2)+n) matrix

    n = sqrt(size(m(1,:),2));
    m_combined = zeros(size(m,1), (((n^2)-n)/2)+n);
    for i=1:size(m,1)
        m_combined(i,:) = combineSymmetricCoefficientsInRow(m(i,:));
    end
end

