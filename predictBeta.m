%% Beta prediction
% This function calculates the one-step-ahead realized beta prediction
% obtained under the model of Moura and Noriller (2019)
%
% Input:
% - k: d.o.f. parameter k
% - sig: covariance scaling parameter
% - B: (q x T) array of realized betas
% 
% Output:
% - Bhat: Beta prediction for T+1

function Bhat = predictBeta(k, sig, B)

[q, T] = size(B);
Q = 1/sig;
lam = k/(k+1);

Iq = eye(q);
% initialization t = 1
b_tt = B(:,1);
N_tt = 1;
S_tt = Iq;

% filtering t = 2:T
for tt = 2:T
    b_tl = b_tt;
    N_tl = 1/(1/(lam*N_tt) + 1/Q);
    S_tl = S_tt*lam*(k+1)/k;
    
    e_t = B(:,tt) - b_tl;
    N_tt = N_tl + 1;
    b_tt = (b_tl*N_tl + B(:,tt))/N_tt;
    S_tt = S_tl*k/(k+1) + e_t*e_t'*(1-1/N_tt)/(k+1);
end

Bhat = b_tt;
end