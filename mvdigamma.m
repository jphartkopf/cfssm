%% Multivariate Digamma (Psi) function
% Digamma_d(x) = derivative of Gamma_d(x) w.r.t x
% with log Gamma_d(x) = d(d-1)/4 log pi + sum_(j=1)^d log Gamma(x+(1-j)/2)
%
% Input:
% - x: argument (scalar real)
% - d: dimension (scalar integer)
% Output:
% - P: value of the logarithmic multivariate Digamma function
%
% Function does minimal input checking, so be careful!

function P = mvdigamma(x, d)
dd = (1:d);
P = sum(psi(x + (1-dd)/2));
end