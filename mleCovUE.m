%% Maximum Likelihood estimation of the UE model
% This function estimates the parameters of the UE model (Windle and
% Carvalho, 2014; Hartkopf, 2020) fitted to the time series of realized 
% covariance matrices under several restrictions.
%
% Input:
% - C: (m x m x T) array of realized covariance matrices
% - par0: vector of initial values for the numerical optimization routine,
% par0 = [n0, k0, lam0]. lam0 is optional for most cases
% - tau1: number of observations set aside for estimation of S0
% - restriction: parameter restriction for the smoothing parameter lambda.
% allowed values are {'R0', 'R1', 'R2', 'R3'}
% 
% Output:
% - n: ML estimate of d.o.f. parameter n
% - k: ML estimate of d.o.f. parameter k
% - lam: ML estimate of smoothing parameter lambda

function [n, k, lam] = mleCovUE(C, par0, tau1, restriction)
[m, ~, T] = size(C); % dimensions of C

opt = optimoptions('fmincon', ...
        'GradObj', 'off', ...
        'Display', 'iter', ...
        'Algorithm', 'sqp', ...
        'maxfuneval', 1e6, ...
        'maxiter', 1e6); 

switch restriction
    case 'R0'
        % no restriction (R0)
        paropt = fmincon(@llh, par0, [], [], [], [], [m+1 m+2 eps], [inf inf 1-eps], [], opt, C, m, T, tau1);
        n = paropt(1);
        k = paropt(2);
        lam = paropt(3);
    case 'R1'
        % EWMA restriction (R1)
        % lam = 1/( 1 + k/(n-m-1) );
        paropt = fmincon(@llhEWMA, par0, [], [], [], [], [m+1 m+2], [], [], opt, C, m, T, tau1);
        n = paropt(1);
        k = paropt(2);
        lam = 1/( 1 + n/(k-m-1) );
    case 'R2'
        % RW restriction (R2)
        % lam = n/(n+k);
        paropt = fmincon(@llhRW, par0, [], [], [], [], [m+1 m+2], [], [], opt, C, m, T, tau1);
        n = paropt(1);
        k = paropt(2);
        lam = k/(k+n);
    case 'R3'
        % psi restriction (R3)
        % lam = exp(E(log(det(U)))/p)
        paropt = fmincon(@llhPSI, par0, [], [], [], [], [m+1 m+2], [inf inf], [], opt, C, m, T, tau1);
        n = paropt(1);
        k = paropt(2);
        lam = exp((-mvdigamma((k+n)/2,m)+mvdigamma(k/2,m))/m);
end

end


%% log likelihood for UE model without anC restriction (R0) on lambda
function f = llh(par, C, m, T, tau1)
n = par(1);
k = par(2);
lam = par(3);

f = (T-tau1)*( mvgammaln((n+k)/2,m) - mvgammaln(n/2,m) - mvgammaln(k/2,m) );

lam_ = permute( lam.^((tau1-1):-1:0) , [1 3 2] );
S = sum(lam_.*C(:,:,1:tau1),3); % S0

for tt = tau1+1:T
    
    S = lam*S;
    f = f + 0.5*(n-m-1)*logdet(C(:,:,tt)) + 0.5*k*logdet(S);
    S = S + C(:,:,tt);
    f = f - 0.5*(k+n)*logdet(S);
    
end

f = - f;

end


%% likelihood for UE model with EWMA (R1) restriction
function f = llhEWMA(par, C, m, T, tau1)
n = par(1);
k = par(2);
lam = 1/( 1+n/(k-m-1) );


f = (T-tau1)*( mvgammaln((n+k)/2,m) - mvgammaln(n/2,m) - mvgammaln(k/2,m) );

lam_ = permute( lam.^((tau1-1):-1:0) , [1 3 2] );
S = sum(lam_.*C(:,:,1:tau1),3); % S0

for tt = tau1+1:T
    
    S = lam*S;
    f = f + 0.5*(n-m-1)*logdet(C(:,:,tt)) + 0.5*k*logdet(S);
    S = S + C(:,:,tt);
    f = f - 0.5*(k+n)*logdet(S);
    
end

f = - f;

end


%% log likelihood for UE model with RW (R2) restriction
function f = llhRW(par, C, m, T, tau1)
n = par(1);
k = par(2);
lam = k/( k+n );


f = (T-tau1)*( mvgammaln((n+k)/2,m) - mvgammaln(n/2,m) - mvgammaln(k/2,m) );

lam_ = permute( lam.^((tau1-1):-1:0) , [1 3 2] );
S = sum(lam_.*C(:,:,1:tau1),3); % S0

for tt = tau1+1:T
    
    S = lam*S;
    f = f + 0.5*(n-m-1)*logdet(C(:,:,tt)) + 0.5*k*logdet(S);
    S = S + C(:,:,tt);
    f = f - 0.5*(k+n)*logdet(S);
    
end

f = - f;

end


%% likelihood for UE model with psi (R3) restriction
function f = llhPSI(par, C, m, T, tau1)
n = par(1);
k = par(2);
lam = exp((-mvdigamma((k+n)/2,m)+mvdigamma(k/2,m))/m);

f = (T-tau1)*( mvgammaln((n+k)/2,m) - mvgammaln(n/2,m) - mvgammaln(k/2,m) );

lam_ = permute( lam.^((tau1-1):-1:0) , [1 3 2] );
S = sum(lam_.*C(:,:,1:tau1),3); % S0

for tt = tau1+1:T
    
    S = lam*S;
    f = f + 0.5*(n-m-1)*logdet(C(:,:,tt)) + 0.5*k*logdet(S);
    S = S + C(:,:,tt);
    f = f - 0.5*(k+n)*logdet(S);
    
end

f = - f;

end

%% end of file