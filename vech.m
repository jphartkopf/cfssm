%% Vectorization of the lower half of the matrix A
% Especially attractive for symmetric matrices, e.g. covariance matrices.
% Some (k x k) matrix is transformed in a m = (k*(k+1)/2 x 1) vector by
% stacking the columns of A beginning with the diagonal entry.
%
% Input:
% - A : (k x k) symmetric matrix
%
% Output:
% - d : (m x 1) vector of the lower triangular of matrix A

function d = vech(A)
d = A( tril( true( size(A) ) ) );
end