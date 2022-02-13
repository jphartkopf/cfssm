%% Trace calculation
% Calculates the trace of square matrix A
%
% Input:
% - A: (m x m) matrix
%
% Output:
% - d: trace of A
%
% Function does minimal input checking, so be careful!

function d = tr(A)
d = sum(diag(A));
end