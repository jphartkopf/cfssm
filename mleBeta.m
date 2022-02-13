%% Maximum Likelihood estimation of the UE model
% This function estimates the parameters of the TVP-VAR (Moura and
% Noriller, 2019; Hartkopf, 2020) fitted to the time series of realized 
% betas.
%
% Input:
% - B: (q x T) array of realized betas
% - par0: initial value of sig0 for the numerical optimization routine
% 
% Output:
% - sig: ML estimate of the scaling parameter sig

function sig = mleBeta(B, par0, k)
[q, T] = size(B); % dimensions of B

opt = optimoptions('fmincon', ...
        'GradObj', 'off', ...
        'Display', 'iter', ...
        'Algorithm', 'sqp', ...
        'maxfuneval', 1e6, ...
        'maxiter', 1e6);
    
sig = fmincon(@llh, par0, [], [], [], [], [eps], [], [], opt, B, T, q, k);
end


%%
function f = llh(par, B, T, q, k)
Q = 1/par(1);
lam = k/(k+1);

Iq = eye(q);
% initialization t = 0
f = 0;
x_tt = B(:,1);
N_tt = 1;
S_tt = Iq;

% filtering t = 1:T
for tt = 2:T
    x_tl = x_tt;
    N_tl = 1/(1/(lam*N_tt) + 1/Q);
    S_tl = S_tt*lam*(k+1)/k;
    
    e_t = B(:,tt) - x_tl;
    N_tt = N_tl + 1;
    x_tt = (x_tl*N_tl + B(:,tt))/N_tt;
    S_tt = S_tl*k/(k+1) + e_t*e_t'*(1-1/N_tt)/(k+1);
    Sig_t = ((1-1/N_tt)*Iq)/(k*S_tl);
    
    f = f - 0.5*q*log((k-q+1)*pi) + gammaln((k+1)/2) - gammaln((k-q+1)/2) ...
        + 0.5*logdet((k-q+1)*Sig_t) - 0.5*(k+1)*log(1+e_t'*Sig_t*e_t);
    
end

f = -f;

end

%% end of file
