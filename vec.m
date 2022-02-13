%% Vectorization of matrix A
% Returns vec(A) by stacking all columns of A in a vector
%
% Input:
% - A : (m x k) matrix
% 
% Output:
% - d : (m*k x 1) vector

function d = vec(A)
d = A(:);
end