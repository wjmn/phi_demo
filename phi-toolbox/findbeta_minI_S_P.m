function [absI_s_p] = findbeta_minI_S_P(beta,C_D_beta1_inv,Cov_X_inv,CovXt,S_left,S_right,I_s_p_Const_part,S)

C_D_beta_inv = beta*C_D_beta1_inv;

%     Q
C_AB_inv = Cov_X_inv + C_D_beta_inv;

%%%
C_AB = inv(C_AB_inv);
%     C_AB = pinv(C_AB_inv);

C_AB_p = -C_AB*C_D_beta1_inv*C_AB;

C_Xd_inv_p = S - beta*S_left*(2*C_AB + beta*C_AB_p)*S_right;

I_s_p = 1/2*(-trace(C_AB_inv*C_AB_p) + trace(CovXt*C_Xd_inv_p)) + I_s_p_Const_part;

% fprintf('iter=%d beta=%f I_s_p=%f\n',iter, beta,I_s_p);

absI_s_p = abs(I_s_p);