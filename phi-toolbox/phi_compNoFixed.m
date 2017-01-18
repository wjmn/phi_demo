function [phi_star phi_star_fixed phi_MI phi_MI_fixed MI MI_fixed phi_H H_cond beta betaFixed] = phi_compNoFixed(CovXt,CovXtXtau,CovXtau,CovXtauFixed,Z,beta_init)

%-----------------------------------------------------------------------
% FUNCTION: phi_comp.m
% PURPOSE:  calculate practical measures of phi given covariances of data
% 
% INPUTS:   
%           Cov_Xt: covariance of data at time t (present state)
%           Cov_XtXtau: cross-covariance of X(t) (present state) and X(t-tau) (past state)
%           Cov_Xtau: covariance of data at time t-tau (past)
%           Cov_fixed: the covariance of the distribution that is to be used (forced) for the past. Used for
%           phi_star_fixed, phi_MI_fixed and I_fixed
%           Z: partition of each channel (default: atomic partition)
%           beta_init: initial value of beta (default: beta_int=1)
%
% OUTPUT:
%           phi_star: phi based on mismatched decoding
%           phi_star_fixed: phi based on mismatched decoding, when
%           Cov_fixed is reinstated as the covariance for the past
%           phi_MI: Barrett and Seth's measure (2011) of phi based on mutual information (Phi E in paper)
%           phi_MI_fixed: same as phi_MI but when the Cov_fixed is used as
%           the covariance for the past (a little different from Phi_DM in
%           Barrett and Seth 2011
%           phi_H: Barrett and Seth's measure (2011) of phi based on conditional entropy (Phi E^-)
%           beta: numerically obtained beta value
%-----------------------------------------------------------------------
%num of elements/channels
N = size(CovXt,1);

%some preliminary checks

%not enough args
if nargin < 2
    error 'not enough input args';
end

%use stationary case - p(past)=p(present)
if nargin < 3
    CovXtau=CovXt;
end

%whetehr to calculate fixed cov ver
% if nargin < 4
    %the only difference is that we will not calculate the gradient decent
    %and set phi_star_fixed, phi_MI_fixed and I_fixed to NaNs 
    CovXtauFixed=eye(N);
    calc_fixed=false;
% else
%     calc_fixed=true;
% end

% %force atmoic partition
% if nargin < 5
%     Z = 1: 1: N;
% end
% 
% %start gradient decent for phi_star with beta=1
% if nargin < 6
     beta_init = 1;
% end

%%% Now start the real work
%% Covariances

% calculate the covariance of the conditional P(Xtau | Xt) {P(past | present)}
CovXtau_Xt = Cov_cond(CovXtau, CovXtXtau', CovXt); 

% calculate the covariance of the conditional P(Xt | Xtau) {P(present | past)}
CovXt_Xtau = Cov_cond(CovXt, CovXtXtau, CovXtau); 

% for fixed prior we use the conditional P(Xt | Xtau) and reinstate Cov_fixed as the
% prior. For this we will need the new conditional Cov_Xfixed_Xt will have the following
% covariance (this is a standard result for applying bayes rule to gaussians
% (e.g. pg 93 of C M Bishop, 2006, Pattern Recognition and Machine
% Learning, also in Barrett Seth's paper))

% First get 'A' (analguous to the connectivity matrix in the stationary
% case)
A=CovXtXtau*CovXtau^-1;

%now to get the conditional over the past when its fixed, given the present
% P(XtauFixed|Xt)
CovXtauFixed_Xt=(CovXtauFixed^-1+A'*(CovXt_Xtau^-1)*A)^-1;

%we will also need the 'new' covariance of the present for phi_star
CovXtFixed=CovXt_Xtau+A*CovXtauFixed*A';

%Importantly, we now have a different cross covariance Cov(Xt,XtauFixed:
CovXtXtauFixed=A*CovXtauFixed;

%% Entropies in the whole system

%first, conditional entropy in the whole system
H_cond = H_gauss(CovXtau_Xt,-1);

% if you are interested in checking the other way of calculating MI then
% you can uncomment the line below. this is not used for the calc of any of
% the outputs, just a check

% % % H_cond2 = H_gauss(CovXt_Xtau,-1);

%second, entropy in the whole system when we reinstated CovXtauFixed on the past
H_condFixed = H_gauss(CovXtauFixed_Xt,-1);

%now for some timely checks. If you see any of these then you probably
%cannot interpret your result meaningfuly...

if isinf(H_cond) == 1
    fprintf('Alert: Infinity Entropy\n')
    H_cond = H_gauss(CovXtau_Xt,-2);
end

if isinf(H_condFixed) == 1
    fprintf('Alert: Infinity Entropy when we used the fixed covariance for the past\n')
    H_condFixed = H_gauss(CovXtauFixed_Xt,-2);
end

if isreal(H_cond) == 0
    fprintf('Alert: Complex Entropy\n')
end

if isreal(H_condFixed) == 0
    fprintf('Alert: Complex Entropy when we used the fixed covariance for the past\n')
end

%% Mutual Info for the whole system

MI = H_gauss(CovXtau,-1) - H_cond; % mutual information

%And if you'd like to verify that the other MI identitiy also works... 
% MI2 = H_gauss(CovXt,-1) - H_cond2; 

%MI when a fixed covariance is used for the parts 
MI_fixed = H_gauss(CovXtauFixed,-1) - H_condFixed;

%And to verify 
% MI_fixed2 = H_gauss(CovXtFixed,-1) - H_cond2;

%And yet another verification, this time using the joint prob

% jointCovFixed=[CovXtFixed CovXtXtauFixed;CovXtXtauFixed' CovXtauFixed];
% MI_fixed3=H_gauss(CovXtFixed,-1)+H_gauss(CovXtauFixed,-1)-H_gauss(jointCovFixed,-1);

%% START CALCULATIONS FOR THE PARTS

% number of parts for this partition
N_c = max(Z); 

%to hold the entropies for the parts for the past P(Mt)
H_p = zeros(N_c,1);
%to hold the entropies of the conditional for the parts P(Mtau|Mt) 
H_cond_p = zeros(N_c,1);
%to hold the mutual information  for the parts MI(Mtau,Mt)
MI_p = zeros(N_c,1);

%now for the same but for the fixed past case
%
H_fixed_p = zeros(N_c,1);
%
H_condFixed_p = zeros(N_c,1);
%
MI_fixed_p = zeros(N_c,1);

%grab the different qunatities for the parts, as desribed by the partition
M_cell = cell(N_c,1);
for i=1: N_c
    M_cell{i} = find(Z==i);
end


%% First do the calculations for Barrett and Seth's measure, Phi_MI and Phi_MI fixed, 2011

% get the different covariances and entropies for the parts

for i=1: N_c
    M = M_cell{i};
    
    % covariance of the part in the past
    %Cov(P(Mtau))
    CovXtau_p = CovXtau(M,M);
    
    %covariance of the part in the present
    %Cov(P(Mt))
    CovXt_p = CovXt(M,M);
 
    %and for the fixed version
    
    %Cov(P(MtauFixed))
    CovXtauFixed_p = CovXtauFixed(M,M);   
    
    %Cov(P(MtauFixed))
    CovXtFixed_p = CovXtFixed(M,M); 
   
    %cross covariance present,past for the parts
    CovXtXtau_p = CovXtXtau(M,M);
    
    %same but fixed
    CovXtXtauFixed_p = CovXtXtauFixed(M,M);
    
    %%% the covariance of the conditionals 
    %cond covariance past_present for the parts Cov(P(Mtau|Mt))  
    CovXtau_Xt_p = Cov_cond(CovXtau_p,CovXtXtau_p',CovXt_p);
    
    %and the other one, only neccesary for checking
    %Cov(P(Mt|Mtau))
%     CovXt_Xtau_p = Cov_cond(CovXt_p,CovXtXtau_p,CovXtau_p);
    
    %%% the covariance of the conditionals when fixed
    CovXtauFixed_XtFixed_p = Cov_cond(CovXtauFixed_p,CovXtXtauFixed_p',CovXtFixed_p);
    
    %Cov(P(Mt|Mtau))
    %only neccessary for checking
%     CovXtFixed_XtauFixed_p = Cov_cond(CovXtFixed_p,CovXtXtauFixed_p,CovXtauFixed_p);
          
    %% Entropy calc
    
    % cond entropy of  P(Mtau|Mt)
    H_cond_p(i) = H_gauss(CovXtau_Xt_p);
    
    % entropy of past for parts
    H_p(i) = H_gauss(CovXtau_p);
    
    %mutual information measured
    MI_p(i) = H_p(i) - H_cond_p(i);
    
    %Of coures the result above is the same as
%     MI_p2=H_gauss(CovXt_p,-1)-H_gauss(CovXt_Xtau_p);
        
    % now for the fixed past case
    
    % cond past|present when past is fixed
    H_condFixed_p(i) = H_gauss(CovXtauFixed_XtFixed_p);
    
    % ent of past when it is fixed
    H_fixed_p(i) = H_gauss(CovXtauFixed_p);  
    
    %and the mutual info 
    MI_fixed_p(i) = H_fixed_p(i) - H_condFixed_p(i);
    
    %of course we could have also done the above as 
%     MI_fixed_p2 = H_gauss(CovXtFixed_p,-1) - H_gauss(CovXtFixed_XtauFixed_p);
    
    %OR
%     subJointCov=[CovXtFixed_p CovXtXtauFixed_p;CovXtXtauFixed_p' CovXtauFixed_p];
%     MI_fixed_p3 = H_gauss(CovXtFixed_p,-1) + H_gauss(CovXtauFixed_p,-1) - H_gauss(subJointCov);
     
    
      
      
    
end

%% fincal calc of phi_MI (and stochastic interaction)

%normalisation 
K = (N_c-1)*min(H_p); % normalization

%normalisation for fixed case 
K_fixed = (N_c-1)*min(H_fixed_p); % normalization

%% 
phi_H(1) = sum(H_cond_p) - H_cond; % entropy-based effective information
phi_H(2) = phi_H(1)/K;

phi_MI(1) = MI - sum(MI_p); % mutual information-based effective information
phi_MI(2) = phi_MI(1)/K;

if calc_fixed
    
    phi_MI_fixed(1) = MI_fixed - sum(MI_fixed_p); % mutual information-based effective information
    phi_MI_fixed(2) = phi_MI_fixed(1)/K_fixed;

else
% if we weren't interested in fixed case then just make it nan    
    phi_MI_fixed(1) = nan;
    phi_MI_fixed(2) = nan;
end

%% Calculation of phi_atar 

%the following hold the block diagnoal matricies required for the phi_star calculation
CovXtau_D = zeros(N,N);
CovXt_D = zeros(N,N);
CovXtXtau_D = zeros(N,N);
C_D_cond = zeros(N,N);
H_cond_D = 0;


%and for the fixed case
CovXtauFixed_D = zeros(N,N);
CovXtFixed_D = zeros(N,N);
CovXtXtauFixed_D = zeros(N,N);
C_D_condFixed = zeros(N,N);
H_condFixed_D = 0;

%these are the same for the fixed and non fixed case


%and we need these again for the parts
for i=1: N_c
    
     M = M_cell{i};
     
     % covariance of the part in the past
    %Cov(P(Mtau))
    CovXtau_p = CovXtau(M,M);
    
    %covariance of the part in the present
    %Cov(P(Mt))
    CovXt_p = CovXt(M,M);
 
    %and for the fidxed fixed version
    
    %Cov(P(MtauFixed))
    CovXtauFixed_p = CovXtauFixed(M,M);   
    
    %Cov(P(MtauFixed))
    CovXtFixed_p = CovXtFixed(M,M); 
   
    %cross covariance present,past for the parts
    CovXtXtau_p = CovXtXtau(M,M);
    
    %same but fixed
    CovXtXtauFixed_p = CovXtXtauFixed(M,M);
    
    %%% the covariance of the conditionals 
    %cond covariance past_present for the parts Cov(P(Mtau|Mt)), only needed for checks  
%     CovXtau_Xt_p = Cov_cond(CovXtau_p,CovXtXtau_p',CovXt_p);
    
    %and the other one, only neccesary for checking
    %Cov(P(Mt|Mtau))
    CovXt_Xtau_p = Cov_cond(CovXt_p,CovXtXtau_p,CovXtau_p);
    
    %%% the covariance of the conditionals when fixed, only needed for
    %%% checks
%     CovXtauFixed_XtFixed_p = Cov_cond(CovXtauFixed_p,CovXtXtauFixed_p',CovXtFixed_p);
    
    %Cov(P(Mt|Mtau))

    %Cov(P(Mt|Mtau))
    CovXtFixed_XtauFixed_p = Cov_cond(CovXtFixed_p,CovXtXtauFixed_p,CovXtauFixed_p);    
    
    %%%the rest are for the calculation of phi_star    
     
    H_cond_D = H_cond_D + H_gauss(CovXt_Xtau_p);
    CovXtau_D(M,M) = CovXtau_p;
    CovXtXtau_D(M,M) = CovXtXtau_p;
    CovXt_D(M,M) = CovXt_p;
    C_D_cond(M,M) = CovXt_Xtau_p;
    
    
    %fixed past stuff
    H_condFixed_D=H_condFixed_D + H_gauss(CovXtFixed_XtauFixed_p);
    CovXtauFixed_D(M,M) = CovXtauFixed_p; %
    CovXtFixed_D(M,M) = CovXtFixed_p;
    CovXtXtauFixed_D(M,M) = CovXtXtauFixed_p;
    C_D_condFixed(M,M) = CovXtFixed_XtauFixed_p;
    
    
end

%% final calcs of phi_star

%and calculate I*
%[MI_star beta iter] = calculateIStarByGradientAscent(C_D_cond, CovXt, CovXtau, CovXtau_D, CovXtXtau_D, H_cond_D, beta_init);
[MI_star beta iter] = calculateIStarByFminsearch(C_D_cond, CovXt, CovXtau, CovXtau_D, CovXtXtau_D, H_cond_D, beta_init);

%for the fixed case
if calc_fixed
    [MI_fixed_star betaFixed iter] = calculateIStarByGradientAscent(C_D_condFixed, CovXtFixed, CovXtauFixed, CovXtauFixed_D,  CovXtXtauFixed_D, H_condFixed_D, beta_init);
else
    MI_fixed_star=nan; beta=nan; iter=nan;
end

K = (N_c-1)*min(H_p); % normalization

%and for the fixed case 
K_fixed = (N_c-1)*min(H_fixed_p); % normalization

%% 

phi_star(1) = MI - MI_star; % our proposed effective information
phi_star(2)  = phi_star(1)/K;

if calc_fixed
    
    phi_star_fixed(1) = MI_fixed - MI_fixed_star; % our proposed effective information
    phi_star_fixed(2)  = phi_star_fixed(1)/K_fixed;

else
    
    phi_star_fixed(1) = nan;
    phi_star_fixed(2)  = nan; 
end

%fprintf('beta=%f I_s=%f\n',beta,I_s);
%fprintf('I=%f phi_star=%f phi_H=%f phi_MI=%f\n',I,phi_star(1),phi_H(1),phi_MI(1));

end
