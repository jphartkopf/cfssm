%% Logarithmic determinant
% This function calculates the logarithmic determinant of some pos.def.
% square matrix A via the Cholesky decomposition of A. In some cases if the
% Cholesky of A cannot be computed due to numerical irregularities, the log
% of the determinant is calculated instead.
%
% Input:
% - A: (m x m) pos. def. matrix
%
% Output:
% - d: logarithmic determinant of A
%
% Function does minimal input checking, so be careful!

function d = logdet(A)
try
    C = chol(A);
    d = 2*sum(log(diag(C)));
catch
    d = log(det(A));
end