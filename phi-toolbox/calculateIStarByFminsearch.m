function [I_s beta iter] = calculateIStarByFminsearch(C_D_cond, CovXt, CovXtau, CovXtau_D, CovXtXtau_D, H_cond_D, beta_init)

% Gradient ascent for beta for phiStar calc

N = size(CovXtau,1);

%%%
%inv of CovD_present_past 
S = inv(C_D_cond);

%inv of CovX
Cov_X_inv = inv(CovXtau);

%second term of Q
C_D_beta1_inv = CovXtau_D\CovXtXtau_D'*S*CovXtXtau_D/CovXtau_D;


%parts of the second term of R
S_left = S'*CovXtXtau_D/CovXtau_D;
S_right = CovXtau_D\CovXtXtau_D'*S;
I_s_p_Const_part = 1/2*(logdet(C_D_cond) + N*log(2*pi)) - H_cond_D;

[beta,~,~,output] = fminsearch(@(beta_init) findbeta_minI_S_P(beta_init,...
                                                 C_D_beta1_inv,...
                                                 Cov_X_inv,...
                                                 CovXt,...
                                                 S_left,...
                                                 S_right,...
                                                 I_s_p_Const_part,...
                                                 S),beta_init);
iter = output.iterations;


C_D_beta_inv = beta*C_D_beta1_inv;

%     Q
C_AB_inv = Cov_X_inv + C_D_beta_inv;

%%%
C_AB = inv(C_AB_inv);
%     C_AB = pinv(C_AB_inv);

C_AB_p = -C_AB*C_D_beta1_inv*C_AB;

%second term of eqn 29 in paper
norm_t = -1/2*logdet(C_AB) + 1/2*logdet(CovXtau) + beta/2*logdet(C_D_cond) ...
    + N*beta/2*log(2*pi);
%first term of eqn 29
C_Xd_inv = beta*S - beta^2*S_left*C_AB*S_right;
trace_t = 1/2*trace(CovXt*C_Xd_inv);

%trace_t is first term on paper, norm_t is second (sign changed), beta*H_cond_D is last term of eqn 29
I_s = norm_t + trace_t - beta*H_cond_D;

end

