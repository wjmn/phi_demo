function detX = logdet(X,Const)
%% compute log of determinant of X

if nargin < 2
    Const = 0;
end

X = X/10^Const;

N = size(X,1);
detX = log(det(X)) + Const*N*log(10);