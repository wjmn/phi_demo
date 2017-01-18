function [Cov_X Cov_XY Cov_YX Cov_Y] = Cov_comp_shrink(X,tau,op_detrend,time_w)

%-----------------------------------------------------------------------
% FUNCTION: Cov_comp.m
% PURPOSE:  estimate the covariance matrix of data by shrinckage estimation
% (Schafer & Strimmer, 2005)
% 
% INPUTS:   
%           X: time series data (in form channels x samples)
%           tau: time difference between past and present
%           op_detrend: detrend data (=1) or not (=0)
%           time_w: length of time bin over which data is detrended
% OUTPUT:
%           Cov_X: covariance matrix of X
%           Cov_XY: cross-covariance of X(t-tau) (past state) and X(t) (present state)
%
% last updated 02/10/2014
%-----------------------------------------------------------------------

N = size(X,1);
T = size(X,2);
if nargin < 3
    op_detrend = 0;
end
if nargin < 4
    time_w = 1000;
end

if nargin > 1
    if op_detrend == 0
        t_range1 = 1:1: T-tau;
        t_range2 = tau+1: 1: T;
        Z = [X(:,t_range2); X(:,t_range1)];
        [S_star lambda] = covshrinkKPM(Z');
        %fprintf('lambda=%f\n',lambda);
        
        Cov_X = S_star(1:N,1:N);
        Cov_XY = S_star(1:N,N+1:2*N);
        Cov_YX = S_star(N+1:2*N,1:N);
        Cov_Y = S_star(N+1:2*N,N+1:2*N);
    else
        t_max = floor(T/time_w);
        Cov_X = zeros(N,N);
        Cov_XY = zeros(N,N);
        Cov_Y = zeros(N,N);
        for t=1: t_max
            st = 1+(t-1)*time_w;
            en = t*time_w;
            t_range = st:en;
            X_p = zeros(N,time_w);
            for i=1: N
                X_p(i,:) = cca_detrend(X(i,t_range)); % detrend
            end
            
            t_range1 = 1: 1: time_w-tau;
            t_range2 = tau+1: 1: time_w;
            Z = [X_p(:,t_range2); X_p(:,t_range1)];
            [S_star lambda] = covshrinkKPM(Z');
            %fprintf('lambda=%f\n',lambda);
            Cov_X = Cov_X + S_star(1:N,1:N)/t_max;
            Cov_XY = Cov_XY + S_star(1:N,N+1:2*N)/t_max; % time-lagged covariance matrix
            Cov_Y = Cov_Y + S_star(N+1:2*N,N+1:2*N)/t_max;
        end
        Cov_YX = Cov_XY';
    end
else
    Cov_X = Cov_comp_sh(X);
end

end