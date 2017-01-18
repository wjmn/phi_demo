function [phi_star_vec phi_star_fixed_vec phi_MI_vec phi_MI_fixed_vec MI_vec MI_fixed_vec H_vec] = phi_time(X,tau,movingwin,Z,options, s_filename,covCalcType)

%-----------------------------------------------------------------------
% FUNCTION: phi_time.m
% PURPOSE:  calculate the time course of phi given time series data
%   
% 
% INPUTS:   
%           X: time series data (in form channels x samples)
%           tau: time difference between past and present (ms)
%           movingwin: in the form [T step] where "T" is the length of time
%                                   window and "step" is the length of step size (ms)
%           Z: partition with which phi is computed
%           s_filename: filename under which results are saved
%           options(1): save results, 0 (do not save) or 1 (save)
%           options(2): parallel computing, 0 (not parallel), 1 (parallel
%           computing with the maximum cores), n (>1) (parallel computing with n cores)
% 
% OUTPUT:
%           phi_star_vec: time course of phi_star, normalised and unnormalised
%           phi_star_fixed_vec: time course of phi_star when a high entropy gaussian is used as teh distrubtion over the past
%           phi_MI_vec: time course of phi_MI, normalised and unnormalised
%           phi_MI_fixed_vec: time course of phi_MI when a high entropy gaussian is used as teh distrubtion over the past
%           MI_vec: time course of mutual information
%           MI_fixed_vec: time course of mutual information
%`          H_vec: mesaured entropy of the distribution at X^t-tau 
%
%-----------------------------------------------------------------------

op_save = options(1); % save
op_parallel = options(2); % parallel computing

%is cov calc type has not been specified
if nargin < 7
    covCalcType='sample';
end
%% parallel computing
% isOpen = matlabpool('size');
% if isOpen == 0 && op_parallel > 0
%     if op_parallel == 1
%         s  = ['matlabpool'];
%     else
%         s = ['matlabpool ' int2str(op_parallel)];
%     end
%     eval(s);
% end
% 
% if isOpen ~= 0 && op_parallel == 0
%     matlabpool close;
% end

%% split the data into windows
T = movingwin(1); % the length of time window
step = movingwin(2); % the length of step size

st_length = floor((size(X,2)-T)/step)+1;
st_vec = 0: step: (st_length)*step; % time vector


X_cell = cell(st_length,1);
for st_i=1: st_length
    st = st_vec(st_i);
    t_range = st+1:1:st+T; % time rage over which phi is computed
    X_cell{st_i} = X(:,t_range);
end
% clear X;

%hold results
phi_star_vec = zeros(st_length,2); % time course of phi_star
phi_star_fixed_vec = zeros(st_length,2); % 
phi_MI_vec = zeros(st_length,2); % time course of phi_MI
phi_MI_fixed_vec = zeros(st_length,2); % 
MI_vec = zeros(st_length,1); % % time course of mutual infomration
MI_fixed_vec = zeros(st_length,1); % 
H_vec = zeros(st_length,1); % time course of entropy
% H_cond_vec = zeros(st_length,1); % time course of conditional entropy

%% fixed covariance
N=size(X_cell{1},1); %num channels
%use identitiy covariance as the covariance over the past
CovXtauFixed=100*eye(N);

%for gradient ascent
beta_init=1;


%% calculate integrated information
for st_i=1: st_length
    
    if mean(mean(corr(X_cell{st_i}')))>0.95
       disp('the data is highly correlated... could get complex entropy')
    end
    
    
    switch covCalcType
        
        case 'sample'
            [CovXt CovXtXtau CovXtau] = Cov_comp_sample(X_cell{st_i},tau);
            
        case 'shrink'
            
            [CovXt CovXtXtau CovXtauXt CovXtau] = Cov_comp_shrink(X_cell{st_i},tau);
            
        case 'default'
            
            error 'unknwon cov calc type'
            
    end
    
    %grab entropy, if this is negative then expect bad things to follow
    H_vec(st_i)=H_gauss(CovXtau,-1);
   
    
%     want this for stationary case?
%     [EI_star I H_cond EI_MI EI_H beta] = phi_comp_old(Cov_X,Cov_XY,Z); may

   [phi_star phi_star_fixed phi_MI phi_MI_fixed MI MI_fixed phi_H beta betaFixed] = phi_comp(CovXt,CovXtXtau,CovXtau,CovXtauFixed,Z,beta_init);
%     disp(['iters performed: ' num2str(st_length-st_i)]);
    
    phi_star_vec(st_i,:) = phi_star;
    phi_star_fixed_vec(st_i,:) = phi_star_fixed;  
    phi_MI_vec(st_i,:) = phi_MI;
    phi_MI_fixed_vec(st_i,:) = phi_MI_fixed;
    MI_vec(st_i) = MI;
    MI_fixed_vec(st_i) = MI_fixed;
%     H_cond_vec(st_i) = H_cond;
end

clear X_cell t_range;

%% save data
if op_save == 1
    s = ['save ' s_filename];
    eval(s);
    fprintf('\n The results are saved at ');
    fprintf(s_filename)
    fprintf('\n')
end

end