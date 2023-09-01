%% This is used to inverly estimate the number of particles given detection rate (size) and number of detection

% 2023/08/19 v1 Created by Zijiang Yang
% 2023/08/20 v2 Posterior re-sampling of nt
% 2023/08/24 v4 Posterior re-sampling of dp

clc, clear, close all

tic

%% Input data
% dp_all = 10:0.1:141.4;   % particle size
dp_all = [50 50];
n_sim = 1000;        % number of simulations
k = 2;               % detected number of particles
n_range = 0:10000;    % actual number of particles

% Exponetial fit for dp distribution
a = 912.9;
b = -0.04293;

%% Create the distribution for dp
dp_prob = a * exp(b * dp_all);
dp_prob = dp_prob / sum(dp_prob); % normalize

%% Randomly sample from the dp distribution
random_samples_dp = randsample(dp_all, n_sim, true, dp_prob);

%% For each random dp sample, get its corresponding posterior and sample from it
random_samples_nt = zeros(1, n_sim);

for i = 1:n_sim
    dp = random_samples_dp(i);
    de = dp / 100;
    p = R_t(de);
    
    % Calculate posterior for the sampled dp
    posterior = zeros(1, length(n_range));
    for j = 1:length(n_range)
        n = n_range(j);
        if n >= k
            posterior(j) = nchoosek(n, k) * p^k * (1-p)^(n-k);
        else
            posterior(j) = 0;
        end
    end
    
    posterior = posterior / sum(posterior);
    
    % Sample from the posterior
    random_samples_nt(i) = randsample(n_range, 1, true, posterior);
end

%% Plot histograms with mean and std
figure;

% Statistics for random_samples_dp
mean_dp = mean(random_samples_dp);
std_dp = std(random_samples_dp);

% Histogram of random_samples_dp
subplot(1, 2, 1);
histogram(random_samples_dp, 'FaceColor', 'b', 'BinWidth', 1);
xlabel('dp');
ylabel('Frequency');
title('Histogram of random samples dp');
hold on; % Keep the histogram on the plot
ylim = get(gca, 'YLim'); % Get the current y-axis limits
plot([mean_dp mean_dp], ylim, '--k', 'LineWidth', 1); % Add vertical line for mean
% Add text for mean and std
text(mean_dp+3, 0.9*ylim(2), sprintf('Mean = %.2f', mean_dp));
text(mean_dp+3, 0.85*ylim(2), sprintf('Std = %.2f', std_dp));
hold off;

% Statistics for random_samples_nt
mean_nt = mean(random_samples_nt);
std_nt = std(random_samples_nt);

% Histogram of random_samples_nt
subplot(1, 2, 2);
histogram(random_samples_nt, 'FaceColor', 'r', 'BinWidth', 10);
xlabel('nt');
ylabel('Frequency');
title('Histogram of random samples nt');
hold on; % Keep the histogram on the plot
ylim = get(gca, 'YLim'); % Get the current y-axis limits
plot([mean_nt mean_nt], ylim, '--k', 'LineWidth', 1); % Add vertical line for mean
% Add text for mean and std
text(mean_nt+3, 0.9*ylim(2), sprintf('Mean = %.2f', mean_nt));
text(mean_nt+3, 0.85*ylim(2), sprintf('Std = %.2f', std_nt));
hold off;

%% Diagnosis statistics
% 1. Calculate summary statistics of random_samples_nt

mean_nt = mean(random_samples_nt);
std_nt = std(random_samples_nt);
median_nt = median(random_samples_nt);

% Using KDE to estimate the density of random_samples_nt
[estimated_density, sample_points] = ksdensity(random_samples_nt);

% Finding the mode using the KDE
[~, index] = max(estimated_density);
mode_kde = sample_points(index);

disp(['Mean of n_t: ', num2str(mean_nt)]);
disp(['Standard deviation of n_t: ', num2str(std_nt)]);
disp(['Median of n_t: ', num2str(median_nt)]);
disp(['Mode of n_t using KDE: ', num2str(mode_kde)]);

% 2. Compute cumulative statistics
n_samples = length(random_samples_nt);

cumulative_means = cumsum(random_samples_nt) ./ (1:n_samples);
cumulative_stds = arrayfun(@(x) std(random_samples_nt(1:x)), 1:n_samples);
cumulative_medians = arrayfun(@(x) median(random_samples_nt(1:x)), 1:n_samples);

% For cumulative mode using KDE, it's computationally intensive, so this is just an example for the first 1000 samples
cumulative_modes_kde = arrayfun(@(x) kde_mode(random_samples_nt(1:x)), 1:n_samples); % You can adjust this range as needed


% 3. Visualization of cumulative statistics
figure;
% Subplot for cumulative means
subplot(2,2,1);
plot(cumulative_means, 'LineWidth', 1.5);
title('Cumulative Mean of n_t');
xlabel('Sample Size');
ylabel('Mean');

% Subplot for cumulative standard deviations
subplot(2,2,2);
plot(cumulative_stds, 'LineWidth', 1.5);
title('Cumulative Standard Deviation of n_t');
xlabel('Sample Size');
ylabel('Standard Deviation');

% Subplot for cumulative medians
subplot(2,2,3);
plot(cumulative_medians, 'LineWidth', 1.5);
title('Cumulative Median of n_t');
xlabel('Sample Size');
ylabel('Median');

% Subplot for cumulative modes using KDE
subplot(2,2,4);
plot(cumulative_modes_kde, 'LineWidth', 1.5);
title('Cumulative Mode of n_t (using KDE)');
xlabel('Sample Size');
ylabel('Mode');


%% Diagnosis statistics - for 95% CI
% Calculate the upper and lower bounds of the 95% CI using quantiles:
CI_lw = quantile(random_samples_nt, 0.025);
CI_up = quantile(random_samples_nt, 0.975);

% Compute the cumulative CI bounds:
cumulative_CI_lw = arrayfun(@(x) quantile(random_samples_nt(1:x), 0.025), 1:n_samples);
cumulative_CI_up = arrayfun(@(x) quantile(random_samples_nt(1:x), 0.975), 1:n_samples);

% Plot the cumulative CI bounds together with the previously calculated cumulative mode and mean:
figure;

% Subplot for cumulative 95% CI, mode using KDE, and mean
hold on;

plot(cumulative_CI_lw, 'Color', [0 0.8 0.6], 'LineWidth', 1.5); % Dark Cyan
plot(cumulative_CI_up, 'Color', [0 0.8 0.6], 'LineWidth', 1.5); % Dark Cyan
plot(cumulative_means, 'Color', [0.4 0.6 1], 'LineWidth', 1.5); % Light Blue
plot(cumulative_modes_kde, 'Color', [1 0.49 0.5], 'LineWidth', 1.5); % Light Red

legend('95% CI Lower', '95% CI Upper', 'Mean', 'Mode (KDE)');
title('Cumulative 95% CI, Mean, and Mode of n_t');
xlabel('Sample Size');

hold off;


%% Store the data
% To create a table storing the cumulative values of CI_lw, CI_up, means, and modes_kde, you can do:
% cumulative_table = table(cumulative_CI_lw', cumulative_CI_up', cumulative_means', cumulative_modes_kde', ...
%     'VariableNames', {'CI_lw', 'CI_up', 'means', 'modes_kde'});

cumulative_table = table(cumulative_means', cumulative_stds', cumulative_medians', cumulative_modes_kde', cumulative_CI_lw', cumulative_CI_up', ...
    'VariableNames', {'Mean', 'StandardDeviation', 'Median', 'Mode_KDE', 'CI_Lower', 'CI_Upper'});


overall_mean = mean(random_samples_nt);
overall_std = std(random_samples_nt);
overall_median = median(random_samples_nt);
overall_mode_kde = kde_mode(random_samples_nt);
overall_CI_lw = quantile(random_samples_nt, 0.025);
overall_CI_up = quantile(random_samples_nt, 0.975);

summary_table = table(overall_mean, overall_std, overall_median, overall_mode_kde, overall_CI_lw, overall_CI_up, ...
    'VariableNames', {'Mean', 'StandardDeviation', 'Median', 'Mode_KDE', 'CI_Lower', 'CI_Upper'});


current_time = datestr(now, 'mmddyy_HHMM');
filename = sprintf('R%s_k%d_nSim%d.csv', current_time, k, n_sim);
writetable(cumulative_table, filename);




toc



function p = R_t(de)

    % Check the value of de and compute the area accordingly
    if de > 0 && de <= 1
        p = pi * (de / 2)^2;
    elseif de > 1 && de < sqrt(2)
        p = (pi - 4 * acos(1/de)) * (de / 2)^2 + de * sin(acos(1/de));
    elseif de > sqrt(2)
        p = 1;
        msgbox('検出率 p > 1なのだ。気を付けてにゃ(*╹▽╹*)', 'Warning', 'warn'); % Pop-up window for warning
    else
        p = NaN; % Set to NaN for invalid
        msgbox('検出率 p < 0なのだ。どこかまちがうんだろうにゃ(*╹▽╹*)', 'Error', 'error'); % Pop-up window for error
    end

end


% Helper function for the cumulative KDE mode computation
function mode_value = kde_mode(data)
    [density, points] = ksdensity(data, 'NumPoints', 1000);
    [~, idx] = max(density);
    mode_value = points(idx);
end


