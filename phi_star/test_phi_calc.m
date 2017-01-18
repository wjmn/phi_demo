function phi_flat = test_phi_calc(data_A, data_B, trials_per_bin, T, tau)
%TEST_PHI_CALC TEST ONLY for question regarding negative phi warnings. 
%   This is a TEMPORARY function only for use to test negative phi
%   warnings. 
% data_A and data_B: trial by sample raw ERP
% trials_per_bin: number of trials per bin over which to average covariance
% matrices
% T: time window in which to calculate covariance matrices
% tau: lag between past and current state. 
% RETURNS: phi, a flat list of phi values. The only purpose of this
% function is for determining negative phi warnings and as such the phi
% values are presented flatly (i.e. not categorised by time bin centre).

if size(data_A) ~= size(data_B)
    error('Data A and B are not the same size!');
end

[num_trials, num_samples] = size(data_A); % should be same as data_B

num_trial_bins = floor(num_trials / trials_per_bin);
num_sample_bins = floor(num_samples / floor(T / 2)) - 1;

phi_flat = zeros((num_trial_bins * num_sample_bins), 1);

% For reproducibility, reset the RNG every time this function runs
rng(0, 'twister');

% Shuffle
shuffle_indices = randperm(num_trials);
data_A = data_A(shuffle_indices, :);
data_B = data_B(shuffle_indices, :);

for idx_trial_bin = 1:num_trial_bins
    
    start_trial = (idx_trial_bin - 1) * trials_per_bin + 1;
    
    for idx_sample_bin = 1:num_sample_bins
        
        start_sample = (idx_sample_bin - 1) * floor(T / 2) + 1;

        data_selection_A = data_A(start_trial:(start_trial + trials_per_bin - 1), ...
                                               start_sample:(start_sample + T));
        data_selection_B = data_B(start_trial:(start_trial + trials_per_bin - 1), ...
                                               start_sample:(start_sample + T));

        Cov_X_cum = zeros(2);
        Cov_XY_cum = zeros(2);

        for idx_trial = 1:trials_per_bin

            data_AB = [data_selection_A(idx_trial, :); data_selection_B(idx_trial, :)];

            [Cov_X_curr, Cov_XY_curr] = Cov_comp_shrink(data_AB, tau);

            Cov_X_cum = Cov_X_cum + Cov_X_curr;
            Cov_XY_cum = Cov_XY_cum + Cov_XY_curr;

        end

        Cov_X = Cov_X_cum / trials_per_bin;
        Cov_XY = Cov_XY_cum / trials_per_bin;

        phi = phi_comp(Cov_X, Cov_XY);
        
        % Use non-normalised value for now
        phi_flat((idx_trial_bin - 1) * num_sample_bins + idx_sample_bin, 1) = phi(1);
        
    end
    
end

