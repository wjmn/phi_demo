function [CovXt CovXtXtau CovXtau] = Cov_comp_sample(X,tau)

%-----------------------------------------------------------------------
% FUNCTION: Cov_comp.m
% PURPOSE:  calculate the covariance matrix of data
% 
% INPUTS:   
%           X: time series data (in form channels x samples)
%           tau: time difference between past and present
% 
% OUTPUT:
%           CovXt: covariance matrix of Xt (present)
%           CovXtau: covariance matrix of Xtau (past)
%           CovXtXtau: cross-covariance of X(t-tau) (past state) and X(t) (present state)
%
%-----------------------------------------------------------------------

N = size(X,1);
T = size(X,2);

mstd = 0;
for i=1: N
    mstd = mstd + std(X(i,:));
end
mstd = mstd/N;
X = X/mstd;

M_X = mean(X,2);
for i=1: N
    X(i,:) = X(i,:) - M_X(i);
end

t_range1 = 1:1: T-tau;
t_range2 = tau+1: 1: T;

 
CovXtXtau = X(:,t_range1)*X(:,t_range2)'/(T-tau);

% CovXtXtau = CovXtXtau'; 


CovXt=cov(X(:,tau:end)');
CovXtau=cov(X(:,1:end-tau)');

%can have a look at the condition number of these to see if they are almost
%signular, but a negative deterimenant is an even clearer

end