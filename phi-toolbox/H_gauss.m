function H = H_gauss(Cov_X,Const)

%-----------------------------------------------------------------------
% FUNCTION: H_gauss.m
% PURPOSE: calculate entropy under the gaussian assumption
% 
% INPUTS:
%           Cov_X: covariance of data X
%           Const: Cov_X is divided by 10^Const
%
% OUTPUT:
%           H: entropy of X
%-----------------------------------------------------------------------

if nargin < 2
    Const = 0;
end

n = size(Cov_X,1);
H = 1/2*logdet(Cov_X,Const) + 1/2*n*log(2*pi*exp(1));