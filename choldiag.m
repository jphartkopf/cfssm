%% LDL decomposition
% Returns a lower triangular matrix L with unit diagonal elements and a
% diagonal positive matrix D, s.th. L*D*L' = A. 
%
% Input:
% - A: (m x m) pos. def. matrix
%
% Output:
% - L: (m x m) lower triangular matrix with unit diagonal elements
% - D: (m x m) diagonal positive matrix
%
% Function does minimal input checking, so be careful!

function [L, D] = choldiag(A)

L = chol(A,'lower');
D = L.*eye(size(A));
L = L/D;
D = D.^2;

end