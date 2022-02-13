%% Logarithmic multivariate Gamma function
% Gamma_d(x) = pi^(d(d-1)/4) prod_(j=1)^d Gamma(x+(1-j)/2)
% log Gamma_d(x) = d(d-1)/4 log pi + sum_(j=1)^d log Gamma(x+(1-j)/2)
%
% Input:
% - x: argument (scalar real)
% - d: dimension (scalar integer)
% Output:
% - G: value of the logarithmic multivariate Gamma function
%
% Function does minimal input checking, so be careful!

function G = mvgammaln(x, d)
dd = (1:d);
G = d*(d-1)/4*log(pi) + sum(gammaln(x + (1-dd)/2));
end