
%% Data loading
% If not already in path with data:
cd('~/Desktop/asconedemo/');
data_path = '/';    % change to path containing ERP data
addpath(genpath(data_path));
% Please also make sure that the practical_phi toolbox is within your
% data_path as we be using it to calculate phi. 

% File names with ERP data
file_1 = '153_4sessions_li187';
file_2 = '153_4sessions_li294';

% Load full data structures
data_struct_1 = load(file_1);
data_struct_2 = load(file_2);

% Load relevant ERP data into new variables for easy access
data_1 = data_struct_1(1).allERP;
data_2 = data_struct_2(1).allERP;

% Release data structures from namespace
clear('data_struct_1');
clear('data_struct_2');

%% Single-Trial Example

% Covariance matrices------------------------------------------------------

% First specify T and tau parameters
T = 200;
tau = 3;

% Take only from the first set of data
current_data_chn_1 = data_1{1};
current_data_chn_2 = data_2{1};

% Take only from the first trial, samples within time bin
current_bin = [current_data_chn_1(1, 1:T); ...
               current_data_chn_2(1, 1:T)];
           
% Subset bin into values for present state and past state
state_present = current_bin(:, (1 + tau): T);
state_past = current_bin(:, 1:(T - tau));

% Find means and extend dimension 2 to same length as state_present and
% state_past
mean_present = mean(state_present, 2);
mean_past = mean(state_past, 2);

mean_present_array = repmat(mean_present, 1, (T - tau));
mean_past_array = repmat(mean_past, 1, (T - tau));

% Normalisation factor
% Can be defined as either N or N - 1 for N samples. Change the value below
% to set the normalisation factor
normaliser = (T - tau) - 1;

% Find deviances from mean of state_present and state_past
dev_present = state_present - mean_present_array;
dev_past = state_past - mean_past_array;

% Covariance matrix estimation
cov_present = (dev_present * dev_present') / normaliser;
cov_cross = (dev_present * dev_past') / normaliser;

% Phi calculation
% We use phi_gauss from the practical_phi toolbox (Oizumi et al 2016) with
% the assumption that covariance of the past state is the same as
% covariance of the present state. As we are only using two channels, we
% undertake the atomic partition by default.
[phi, I, H] = phi_gauss(cov_present, cov_cross);

% Print values to console to view
% Note that the output phi contains 2 elements; we are only interested in
% the first for this 2-channel subsystem (the second is the normalised phi
% to account for different subsystem sizes).
phi
I
H

%% Phi calculations over all trials
%
% I have saved the single-trial script into a generalisable function named
% calc_single_trial which takes input current_bin and tau. 
% 
% We wish to calculate phi structures for all of the data. We first will
% need to initialise a structure to store all of our phi values. We will
% store our phi values in a 1x5 structure (1 x class of visual stimuli)
% with 4 fields: bin_centres, values_H, values_I and values_phi. Each of 
% values_H,  values_I and values_phi are a trial x time bin array 
% containing the raw phi values for each trial in each time bin. 
% bin_centres is a 1 x time bin array of time bin centres.

processed = struct();


% This rides on the assumption that data_1 and data_2 contain the same and
% corresponding classes.
for class = 1:length(data_1)
    
    current_data_chn_1 = data_1{class};
    current_data_chn_2 = data_2{class};
    
    % bin_centre will be the same for every class since we have the
    % same number of samples. I will, however, include this in every 
    % class - this will be redundant, but I have chosen to do this so all 
    % the needed information is stored within the one structure.
    bin_centres = (T/2):(T/2):(size(current_data_chn_1, 2) - T/2);
    % I have clipped the end by T/2 to ensure we do not iterate beyond the
    % number of samples.
    
    % Initialise matrices to store values of H, I and phi
    num_trials = size(current_data_chn_1, 1);
    num_bins = length(bin_centres);
    values_H = zeros(num_trials, num_bins);
    values_I = zeros(num_trials, num_bins);
    values_phi = zeros(num_trials, num_bins);
    
    % Iterate through every time bin 
    % Take only from the first trial, samples within time bin
    % Nesting for loops like this is generally not recommended. I have done
    % this for convenience but a vectorized approach might be better in the
    % future. 
    for centre = bin_centres
        
        for trial = 1:num_trials
        
            % Find start and end indices for indexing current data
            bin_start = centre - (T/2) + 1;
            bin_end = centre + (T/2);

            current_bin = [current_data_chn_1(trial, bin_start:bin_end); ...
                           current_data_chn_2(trial, bin_start:bin_end)];
            
            % Calculate phi
            [phi_pair, I, H] = calc_single_trial(current_bin, tau);
            
            % Only keep first value of phi (non-normalised)
            phi = phi_pair(1);
            
            % Prints to console if negative phi value obtained
            if phi < 0
                fprintf('Alert: Negative Phi');
            end
            
            % Add to storage matrices
            bin_index = find(bin_centres == centre);
            values_H(trial, bin_index) = H;
            values_I(trial, bin_index) = I;
            values_phi(trial, bin_index) = phi;

        end
    end
    
    % Write the data to the structure
    processed(class).bin_centres = bin_centres;
    processed(class).values_H = values_H;
    processed(class).values_I = values_I;
    processed(class).values_phi = values_phi;
    
end

%% Visualising single-trial phi 
%
% Using imagesc, we visualise a trial x time bin colour map of H, I and
% phi.


fig = figure;
for row = 1:3
    
    % Find minimum and maximum values so colour maps within a row are
    % consistent. Placeholder till I think of a better way to do this.
    H_mins = [];
    H_maxs = [];
    I_mins = [];
    I_maxs = [];
    phi_mins = [];
    phi_maxs = [];
    for i = 1:5
        H_mins = [H_mins min(min(real(processed(i).values_H)))];
        H_maxs = [H_maxs max(max(real(processed(i).values_H)))];
        I_mins = [I_mins min(min(real(processed(i).values_I)))];
        I_maxs = [I_maxs max(max(real(processed(i).values_I)))];
        temp = processed(i).values_phi;
        temp(imag(temp) ~= 0) = NaN;
        phi_mins = [phi_mins min(min((temp)))];
        phi_maxs = [phi_maxs max(max((temp)))];
    end
    
    bounds = [min(H_mins), min(I_mins), min(phi_mins); ...
              max(H_maxs), max(I_maxs), max(phi_maxs)];
    
    for column = 1:5
        
        subplot(3, 5, ((row - 1)* 5 + column));
        
        parameters = {'values_H', 'values_I', 'values_phi'};
        parameter = parameters{row};       
        
        eval(['plot_data = processed(column).' parameter ';']);

        % Remove complex values
        plot_data(imag(plot_data) ~= 0) = NaN;
        
        imagesc(plot_data);
        caxis manual
        caxis([bounds(1, row), bounds(2, row)]);

    end
end

% Save the figure
saveas(fig, 'phi_structures.png');

%% Visualising mean phi

fig_means = figure;
for row = 1:3
    
    for column = 1:5
        
        subplot(3, 5, ((row - 1)* 5 + column));
        
        parameters = {'values_H', 'values_I', 'values_phi'};
        parameter = parameters{row};       
        
        eval(['plot_data = processed(column).' parameter ';']);
        
        % Remove complex values
        plot_data(imag(plot_data) ~= 0) = NaN;
        
        plot_means = nanmean(plot_data, 1);
        plot_stderrs = nanstd(plot_data, 1, 1);

        patch([bin_centres, fliplr(bin_centres)], ...
            [plot_means + plot_stderrs, fliplr(plot_means - plot_stderrs)], ...
            [0.85, 0.85, 0.85], ... % colour
            'EdgeColor', 'None');
        hold on
        plot(bin_centres, plot_means)

    end
end

% Note the out of bounds is because we are using a line graph to display
% discrete time bins. Will think of how to fix later.
saveas(fig_means, 'phi_structure_means.png');

%% TODO
% TODO Add labels to plot above
% TODO plot mean across all trials as a line graph with mean and std
% shading
% TODO possibly bin into 3 trial bins?