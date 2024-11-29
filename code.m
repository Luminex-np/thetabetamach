clc;
clear;

% Define the range of beta values in degrees
beta_deg = linspace(0, 90, 9001);  

% Number of Mach numbers to evaluate
num_Mach = 20;

% Initialize arrays for maximum theta and corresponding beta values
max_theta = zeros(num_Mach, 1);
max_beta = zeros(num_Mach, 1);

% Convert beta values to radians
beta_rad = deg2rad(beta_deg);

% Loop over different Mach numbers
figure;
hold on;
for k = 1:num_Mach
    Mach = 1.1^k; % Mach number
    theta_rad = atan(2 * cot(beta_rad) .* ...
        ((Mach^2 .* sin(beta_rad).^2 - 1) ./ ...
                   (Mach^2 * (1.4 + cos(2 * beta_rad)) + 2)));
    theta_deg = rad2deg(theta_rad); % Convert theta to degrees

    % Plot theta-beta curve
    plot(theta_deg, beta_deg, 'Color', [0.5 0.5 0.5]);  % Gray color for individual curves

    % Find the maximum theta value and corresponding beta value
    [max_theta(k), idx_max] = max(theta_deg);
    max_beta(k) = beta_deg(idx_max);
end

% Plot the curve for maximum theta values
plot(max_theta, max_beta, 'LineWidth', 2, 'DisplayName', ...
    'Max \theta Curve', 'Color', 'b');

% Add labels and grid
xlabel('\theta (degrees)', 'FontSize', 20);
ylabel('\beta (degrees)', 'FontSize', 20);
grid on;
grid minor;
legend show;

% Plot the sonic line
sonic_theta = linspace(0, 45, 500);
sonic_beta = 90 - sonic_theta;
plot(sonic_theta, sonic_beta,'k--','DisplayName', 'Sonic Line','Color','r');

% Set x-axis limit
xlim([0 45]);

% Customize the appearance
title('Theta-Beta Relationship for Various Mach Numbers');
legend('Location', 'best');
hold off;
