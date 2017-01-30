function [phi, I, H] = calc_single_trial(current_bin_1, current_bin_2, tau)
%Calculates phi structures for a single trial of ERP value data
    %
    % INPUTS:
    %   current_bin_1 and 2: a trial x sample array of ERP values for channel
    %   1 and 2
    %   tau: the time lag value to distinguish present and past state
    % OUTPUTS:
    %   phi: pair of non-normalised and normalised phi values
    %   I: mutual information 
    %   H: conditional entropy
    
    % Random number generator seed
    rng(0, 'twister');
    
    [trials_per_bin, ~] = size(current_bin_1);
    
    % Calculating covariance matrices
    Cov_present_cum = zeros(2);
    Cov_cross_cum = zeros(2);
    Cov_past_cum = zeros(2);

    for idx_trial = 1:trials_per_bin

        data_AB = [current_bin_1(idx_trial, :); current_bin_2(idx_trial, :)];

        [Cov_present_curr, Cov_cross_curr, ~, Cov_past_curr] = Cov_comp_shrink(data_AB, tau);

        Cov_present_cum = Cov_present_cum + Cov_present_curr;
        Cov_cross_cum = Cov_cross_cum + Cov_cross_curr;
        Cov_past_cum = Cov_past_cum + Cov_past_curr;

    end

    Cov_present = Cov_present_cum / trials_per_bin;
    Cov_cross = Cov_cross_cum / trials_per_bin;
    Cov_past = Cov_past_cum / trials_per_bin;

    %% Phi calculation
    % We use phi_gauss from the practical_phi toolbox (Haun et al 2016) with
    % the assumption that covariance of the past state is the same as
    % covariance of the present state. As we are only using two channels, we
    % undertake the atomic partition by default.
    [phi, ~, ~, ~, I, ~, ~, H] = phi_comp(Cov_present, Cov_cross, Cov_past);

end