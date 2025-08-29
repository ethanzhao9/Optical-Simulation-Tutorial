%% STED Donut Beam Visualization
% Formula:
%   I(r) = I0 * (r^2 / w^2) * exp(-2*r^2 / w^2)
clear; clc; close all;

% Parameters
I0 = 1;                         % Normalization factor
w_list = [0.5, 1.0, 2.0];       % Different w values (beam waists)

% Grid settings (spatial coordinates)
N = 512;                        % Grid resolution
L = 4;                          % Half-size of field of view (arbitrary units)
x = linspace(-L, L, N);
y = linspace(-L, L, N);
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);          % Radial distance from center

figure('Position', [100 100 1200 600]);

for k = 1:length(w_list)
    w = w_list(k);

    % Donut beam intensity (normalized)
    I = I0 * (R.^2 / w^2) .* exp(-2 * R.^2 / w^2);
    I = I / max(I(:));

    % --- Plot 2D donut intensity map ---
    subplot(2, length(w_list), k);
    imagesc(x, y, I);
    axis image; colormap hot; 
    title(sprintf('2D Donut (w=%.1f)', w));
    xlabel('x (a.u.)'); ylabel('y (a.u.)');

    % --- Plot radial cross-section ---
    % Take profile along central horizontal line (y=0)
    center = round(N/2);
    I_line = I(center,:);
    subplot(2, length(w_list), k + length(w_list));
    plot(x, I_line, 'LineWidth', 1.8);
    grid on;
    xlabel('r (a.u.)'); ylabel('Normalized Intensity');
    title(sprintf('Radial Profile (w=%.1f)', w));
end