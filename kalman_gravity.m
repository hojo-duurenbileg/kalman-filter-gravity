% Kalman filter estimate of gravity from falling-mass position data.
%
% Compares two dynamic models (random constant and random walk) and uses the
% innovation sequence for outlier detection (chi-squared test, 95%).
%
% Author: Khongorzul Duurenbileg, November 2019
%
% Input: two-column text file, time (s) and measured vertical position (m).
close all
clear
clc
data = load('data/sample_falling_mass.txt');
t = data(:, 1);
z_meas = data(:, 2);
% plot
figure(1)
plot(data(:, 1), data(:, 2), '.r');
grid on
title('Measurements of Vertical Position vs. Time');
xlabel('Time (s)');
ylabel('Measured Vertical Position (m)');
Po_ = 0.0001; % gravity variance
Rk = 0.0025; % vertical position variance
xo_ = 9.8;         % gravity
%% Random constant process
PHI = 1;
Qk = 0;
xk_ = PHI * xo_;
Pk_ = PHI * Po_ * PHI + Qk;
x_hat_const = zeros(length(z_meas), 1); % gravity
P_hat_const = zeros(length(z_meas), 1); % gravity variance
vk_const = zeros(length(z_meas), 1); % innovation sequence residual
Cvk_const = zeros(length(z_meas), 1); % innovation sequence VC
T_test_const = zeros(length(z_meas), 1); % global and local test
for i = 1:length(z_meas)
    Hk = 0.5 * t(i)^2;
    Kk = Pk_ * Hk' * inv(Hk * Pk_ * Hk' + Rk); %#ok<MINV>
    xk_plus = xk_ + Kk * (z_meas(i) - Hk * xk_);
    Pk_plus = (1 - Kk * Hk) * Pk_;
    % innovation sequence
    vk_const(i) = z_meas(i) - Hk * xk_;
    Cvk_const(i) = Hk * Pk_ * Hk' + Rk;
    % T-test
    T_test_const(i) = vk_const(i) * (1/Cvk_const(i)) * vk_const(i);
    xk_ = PHI * xk_plus;
    Pk_ = PHI * Pk_plus * PHI' + Qk;
    x_hat_const(i) = xk_;
    P_hat_const(i) = Pk_;
end
outlier_const = T_test_const > 3.84;
% plot
figure(2)
plot(t, x_hat_const, '.b');
grid on
title('Random Constant Process: Gravity vs. Time');
xlabel('Time (s)');
ylabel('Estimated gravity (m/s^{2})');
figure(3)
p = plot(t, vk_const, '-k');
p.LineWidth = 2;
grid on
title('Random Constant Process: Innovation v_{k} vs. Time');
xlabel('Time (s)');
ylabel('Innovation (m)');
%% Random walk process
PHI = 1;
Qk = 0.02 * 0.5;
xk_ = PHI * xo_;
Pk_ = PHI * Po_ * PHI + Qk;
x_hat_walk = zeros(length(z_meas), 1); % gravity
P_hat_walk = zeros(length(z_meas), 1); % gravity variance
vk_walk = zeros(length(z_meas), 1); % innovation sequence residual
Cvk_walk = zeros(length(z_meas), 1); % innovation sequence VC
T_test_walk = zeros(length(z_meas), 1); % global and local test
for i = 1:length(z_meas)
    Hk = 0.5 * t(i)^2;
    Kk = Pk_ * Hk' * inv(Hk * Pk_ * Hk' + Rk); %#ok<MINV>
    xk_plus = xk_ + Kk * (z_meas(i) - Hk * xk_);
    Pk_plus = (1 - Kk * Hk) * Pk_;
    % innovation sequence
    vk_walk(i) = z_meas(i) - Hk * xk_;
    Cvk_walk(i) = Hk * Pk_ * Hk' + Rk;
    % T-test
    T_test_walk(i) = vk_walk(i) * (1/Cvk_walk(i)) * vk_walk(i);
    xk_ = PHI * xk_plus;
    Pk_ = PHI * Pk_plus * PHI' + Qk;
    x_hat_walk(i) = xk_;
    P_hat_walk(i) = Pk_;
end
% plot
figure(4)
plot(t, x_hat_walk, '.b');
grid on
title('Random Walk Process: Gravity vs. Time');
xlabel('Time (s)');
ylabel('Estimated gravity (m/s^{2})');
figure(5)
plot(t, P_hat_walk, '-.k');
grid on
title('Random Walk Process: Pk vs. Time');
xlabel('Time (s)');
ylabel('Pk (m)');
figure(6)
p = plot(t, vk_walk, '-k');
p.LineWidth = 2;
grid on
title('Random Walk Process: Innovation v_{k} vs. Time');
xlabel('Time (s)');
ylabel('Innovation (m)');
figure(7)
p = plot(t, x_hat_const, '-.b');
p.LineWidth = 2;
hold on
p = plot(t, x_hat_walk, '-r');
p.LineWidth = 1;
grid on
title('Kalman Filter Analysis: Gravity vs. Time');
xlabel('Time (s)');
ylabel('Estimated gravity (m/s^{2})');
legend('Random Constant', 'Random Walk')
figure(8)
p = plot(t, vk_const, '-.b');
p.LineWidth = 2;
hold on
p = plot(t, vk_walk, '-r');
p.LineWidth = 1;
grid on
title('Innovation v_{k} vs. Time');
xlabel('Time (s)');
ylabel('Innovation (m)');
legend('Random Constant', 'Random Walk')
