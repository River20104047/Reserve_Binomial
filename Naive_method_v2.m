%% This is used to estimate concentration by Naive approach
% 2023/09/27 By Zijiang Yang
% 2023/09/27 v2: convert output concentration in pc/m3

%% Preparework space
clc, clear, close all

tic

%% Calculations

% Provided Data
V_sea = 850 * 1e3;          % Seawater volume in mL
V_sp = 54;                  % Sample volume in mL
v_sp_i = 0.024255;          % Subsample volume in mL
I = 9;                      % Total Number of Subsamples
n_o = [0 0 1 1 0 0 0 0];    % Observed Particles in each subsample
Q = 10000;                  % Resampling iterations

% Calculations
sum_v_sp = v_sp_i * I;
n_m = sum(n_o);

% Bootstrapping for n_-m
bootstrap_samples = zeros(Q, 1);
for q = 1:Q
    resampled_data = randsample(n_o, length(n_o), true);
    bootstrap_samples(q) = sum(resampled_data);
end
n_hat_m = mean(bootstrap_samples);
quantiles_n_hat_m = quantile(bootstrap_samples, [0.025, 0.975]);

% Concentration calculations
C_m = (V_sp/sum_v_sp) * n_m / V_sea;

C_minus_m_values = (V_sp/sum_v_sp) * bootstrap_samples / V_sea;
quantiles_C_minus_m = quantile(C_minus_m_values, [0.025, 0.975]);

C = (sum_v_sp/V_sp) * C_m + (V_sp - sum_v_sp)/V_sp * mean(C_minus_m_values);
C_values = (sum_v_sp/V_sp) * C_m + (V_sp - sum_v_sp)/V_sp * C_minus_m_values;
quantiles_C = quantile(C_values, [0.025, 0.975]);

% Results
C = C * 1e6; % converting to pieces/m3
quantiles_C = quantiles_C * 1e6; % converting to pieces/m3

C_m = C_m * 1e6; % converting to pieces/m3

C_minus_m = mean(C_minus_m_values) * 1e6; % converting to pieces/m3
quantiles_C_minus_m = quantiles_C_minus_m * 1e6; % converting to pieces/m3

% Display in a Table
T = table({'C'; 'C_m'; 'C_-m'}, [C; C_m; C_minus_m], [quantiles_C(1); NaN; quantiles_C_minus_m(1)], [quantiles_C(2); NaN; quantiles_C_minus_m(2)]);
T.Properties.VariableNames = {'Parameter', 'Estimated Value', '2.5% Quantile', '97.5% Quantile'};

disp(T);


% Calculating the cumulative mean and quantiles for C
cumulative_mean_C = cumsum(C_values)./ (1:Q)';
cumulative_quantiles_C = zeros(Q, 2);

for q = 1:Q
    cumulative_quantiles_C(q, :) = quantile(C_values(1:q), [0.025, 0.975]);
end

% Plotting
figure;
plot(1:Q, cumulative_mean_C * 1e6, 'k', 'LineWidth', 2); % Converting to pieces/m3
hold on;
plot(1:Q, cumulative_quantiles_C(:, 1) * 1e6, '--r'); % Converting to pieces/m3
plot(1:Q, cumulative_quantiles_C(:, 2) * 1e6, '--b'); % Converting to pieces/m3
xlabel('Number of Bootstrap Iterations');
ylabel('Concentration (pieces/m^3)');
title('Cumulative Mean and Quantiles of C');
legend('Cumulative Mean', '2.5% Quantile', '97.5% Quantile');
grid on;
hold off;


toc