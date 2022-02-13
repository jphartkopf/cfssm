%% Covariance prediction
% This function calculates the one-step-ahead realized covariance
% prediction obtained under the UE model of Windle and Carvalho (2014)
%
% Input:
% - n: d.o.f. parameter n
% - k: d.o.f. parameter k
% - lam: smoothing parameter lambda
% - C: (m x m x T) array of realized covariance matrices
% - tau1: number of observations set asied for estimation of S0
% 
% Output:
% - Chat: Covariance prediction for T+1

function Chat = predictCovUE(n, k, lam, C, tau1)

[m,~, T] = size(C); % dimensions of C
lam_ = permute(lam.^((tau1-1):-1:0), [1 3 2]);
S = sum(lam_.*C(:,:,1:tau1),3); % S0

for tt = tau1+1:T
    S = lam*S + C(:,:,tt);
end

Chat = lam*n*S/(k-m-1);

end

%% end of file