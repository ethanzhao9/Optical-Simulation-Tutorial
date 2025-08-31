% Simulate Moiré Fringes to Illustrate a Key Aspect of SIM
% Author: Ethan Zhao

clear; clc; close all;

%% ====== Simulation Parameters ======
grid_size = 512;        % Image grid size (pixels)
[X, Y] = meshgrid(1:grid_size, 1:grid_size); % Coordinate grid

% --- Sample Pattern (Fine Details) ---
% This represents the high-frequency information in your specimen.
% For SIM, this 'freq_sample' would often be a frequency
% that is difficult or impossible for a conventional microscope to resolve.
freq_sample = 22;       % Spatial frequency of the sample (cycles per grid_size)
angle_sample_deg = 0; % Angle of the sample grating lines (degrees from vertical)

% --- Illumination Pattern ---
% This is the structured light pattern projected onto the sample.
% In SIM, 'freq_illumination' is a known, high spatial frequency,
% often chosen to be close to the resolution limit of the microscope.
% The Moiré effect arises from the interaction of freq_sample and freq_illumination.
% The key Moiré frequency will be |freq_sample - freq_illumination|.
freq_illumination = 20; % Spatial frequency of the illumination pattern
angle_illumination_deg = 5; % Angle of the illumination grating lines (degrees)

fprintf('Simulating Moiré fringes...\n');
fprintf('Sample: Frequency = %.1f, Angle = %.1f deg\n', freq_sample, angle_sample_deg);
fprintf('Illumination: Frequency = %.1f, Angle = %.1f deg\n', freq_illumination, angle_illumination_deg);

%% ====== Generate Patterns ======

% Convert angles to radians for trigonometric functions
angle_sample_rad = deg2rad(angle_sample_deg);
angle_illumination_rad = deg2rad(angle_illumination_deg);

% Normalize frequencies to cycles per pixel
k_sample_pix = freq_sample / grid_size; % Spatial frequency in cycles/pixel
k_illum_pix = freq_illumination / grid_size; % Spatial frequency in cycles/pixel

% Create Sample Pattern
% G(x,y) = 0.5 * (1 + cos(2*pi*(k_x*x + k_y*y)))
% k_x = k_pix * cos(angle), k_y = k_pix * sin(angle)
% The term 'X * cos(angle_sample_rad) + Y * sin(angle_sample_rad)'
% effectively rotates the coordinate system for the grating.
sample_pattern = 0.5 * (1 + cos(2 * pi * k_sample_pix * ...
    (X * cos(angle_sample_rad) + Y * sin(angle_sample_rad))));

% Create Illumination Pattern
illumination_pattern = 0.5 * (1 + cos(2 * pi * k_illum_pix * ...
    (X * cos(angle_illumination_rad) + Y * sin(angle_illumination_rad))));

%% ====== Simulate Image Formation (Moiré Effect) ======
% In a simplified model, the observed image is the product of the
% sample's structure and the illumination pattern.
moire_image = sample_pattern .* illumination_pattern;

%% ====== Fourier Transforms ======
% Compute FFTs to see the frequency components.
% The fftshift centers the DC component. We take log(1+abs(FFT)) for visualization.
% Subtracting the mean before FFT helps in reducing the central DC spike's dominance.
compute_fft_display = @(img) fftshift(log(1 + abs(fft2(img - mean(img(:))))));

fft_sample = compute_fft_display(sample_pattern);
fft_illumination = compute_fft_display(illumination_pattern);
fft_moire = compute_fft_display(moire_image);

%% ====== Visualization ======
figure('Position', [100, 100, 1300, 850]);

% --- Spatial Domain ---
h_ax1 = subplot(2, 3, 1);
imagesc(sample_pattern);
colormap(h_ax1, 'parula'); axis image off;
title(sprintf('Sample Pattern\nFreq: %.1f, Angle: %.1f deg', freq_sample, angle_sample_deg));

h_ax2 = subplot(2, 3, 2);
imagesc(illumination_pattern);
colormap(h_ax2, 'parula'); axis image off;
title(sprintf('Illumination Pattern\nFreq: %.1f, Angle: %.1f deg', freq_illumination, angle_illumination_deg));

h_ax3 = subplot(2, 3, 3);
imagesc(moire_image);
colormap(h_ax3, 'parula'); axis image off;
title('Resulting Image (Moiré Fringes)');
xlabel({'Moiré fringes are typically coarser', '(lower frequency)'}, 'FontSize', 10);

% --- Frequency Domain (Fourier Space) ---
% Note on FFT line artifacts: Pure sinusoidal gratings on a discrete grid can
% sometimes produce cross-like line artifacts in their log-scaled FFT display,
% especially if they don't complete an integer number of cycles over the grid.
% The circular markers indicate the true fundamental frequencies.

% Define center of FFT for plotting markers
fft_center = grid_size / 2 + 1;

h_ax4 = subplot(2, 3, 4);
imagesc(fft_sample);
colormap(h_ax4, 'parula'); axis image off;
title('FFT of Sample Pattern');
hold on;
plot(fft_center + freq_sample*cos(angle_sample_rad), fft_center - freq_sample*sin(angle_sample_rad), ...
     'ro', 'MarkerFaceColor','r', 'MarkerSize', 7, 'DisplayName', sprintf('Sample Freq (%.1f)', freq_sample));
plot(fft_center - freq_sample*cos(angle_sample_rad), fft_center + freq_sample*sin(angle_sample_rad), ...
     'ro', 'MarkerFaceColor','r', 'MarkerSize', 7, 'HandleVisibility','off'); % Symmetric component
hold off;
legend('show', 'Location', 'southoutside', 'FontSize', 8);

h_ax5 = subplot(2, 3, 5);
imagesc(fft_illumination);
colormap(h_ax5, 'parula'); axis image off;
title('FFT of Illumination Pattern');
hold on;
plot(fft_center + freq_illumination*cos(angle_illumination_rad), fft_center - freq_illumination*sin(angle_illumination_rad), ...
     'go', 'MarkerFaceColor','g', 'MarkerSize', 7, 'DisplayName', sprintf('Illumination Freq (%.1f)', freq_illumination));
plot(fft_center - freq_illumination*cos(angle_illumination_rad), fft_center + freq_illumination*sin(angle_illumination_rad), ...
     'go', 'MarkerFaceColor','g', 'MarkerSize', 7, 'HandleVisibility','off');
hold off;
legend('show', 'Location', 'southoutside', 'FontSize', 8);

h_ax6 = subplot(2, 3, 6);
imagesc(fft_moire);
colormap(h_ax6, 'parula'); axis image off;
title('FFT of Moiré Image');

% Calculate frequency vector components in cycles per grid (for plotting)
fx_s = freq_sample * cos(angle_sample_rad);
fy_s = freq_sample * sin(angle_sample_rad);
fx_i = freq_illumination * cos(angle_illumination_rad);
fy_i = freq_illumination * sin(angle_illumination_rad);

% Difference frequencies (Moiré)
fx_diff = fx_s - fx_i;
fy_diff = fy_s - fy_i;
freq_diff = sqrt(fx_diff^2 + fy_diff^2);

% Sum frequencies (Moiré)
fx_sum = fx_s + fx_i;
fy_sum = fy_s + fy_i;
freq_sum = sqrt(fx_sum^2 + fy_sum^2);

hold on;
% Original sample peaks (for reference)
plot(fft_center + fx_s, fft_center - fy_s, 'ro', 'MarkerSize', 5, 'DisplayName', sprintf('Orig Sample (%.1f)', freq_sample));
plot(fft_center - fx_s, fft_center + fy_s, 'ro', 'MarkerSize', 5, 'HandleVisibility','off');
% Original illumination peaks (for reference)
plot(fft_center + fx_i, fft_center - fy_i, 'go', 'MarkerSize', 5, 'DisplayName', sprintf('Orig Illum (%.1f)', freq_illumination));
plot(fft_center - fx_i, fft_center + fy_i, 'go', 'MarkerSize', 5, 'HandleVisibility','off');

% Moiré (difference) peaks
plot(fft_center + fx_diff, fft_center - fy_diff, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 9, 'DisplayName', sprintf('Moiré Diff (%.1f)', freq_diff));
plot(fft_center - fx_diff, fft_center + fy_diff, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 9, 'HandleVisibility','off');
% Moiré (sum) peaks
plot(fft_center + fx_sum, fft_center - fy_sum, 'mo', 'MarkerFaceColor','m', 'MarkerSize', 7, 'DisplayName', sprintf('Moiré Sum (%.1f)', freq_sum));
plot(fft_center - fx_sum, fft_center + fy_sum, 'mo', 'MarkerFaceColor','m', 'MarkerSize', 7, 'HandleVisibility','off');
hold off;
legend('show', 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);

sgtitle('Demonstration of Moiré Fringes and their Frequency Components', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('Done.\n');
%% ====== Visualization ======
% Modified to focus on zoomed FFTs and highlight frequency components.

fprintf('Displaying zoomed FFTs of Sample, Illumination, and Moiré Image...\n');

figure('Position', [100, 100, 1500, 550]); % Adjusted for 1x3 layout

% Define center of FFT for plotting markers
fft_center = grid_size / 2 + 1;

% Determine zoom extent for FFT plots
% We want to comfortably see the sample, illumination, and Moiré difference frequencies.
% The sum frequency might be further out.
zoom_radius_pixels = max(freq_sample, freq_illumination) * 1.8; % Zoom factor to ensure visibility

% --- FFT of Sample Pattern ---
h_ax1 = subplot(1, 3, 1);
imagesc(fft_sample);
colormap(h_ax1, 'parula'); axis image off;
title({'FFT of Sample Pattern';'(Zoomed View)'}, 'FontSize', 10);
hold on;
plot(fft_center + freq_sample*cos(angle_sample_rad), fft_center - freq_sample*sin(angle_sample_rad), ...
     'ro', 'MarkerSize', 15, 'LineWidth', 0.5, 'DisplayName', sprintf('Sample Freq (%.1f)', freq_sample));
plot(fft_center - freq_sample*cos(angle_sample_rad), fft_center + freq_sample*sin(angle_sample_rad), ...
     'ro', 'MarkerSize', 15, 'LineWidth', 0.5, 'HandleVisibility','off'); % Symmetric component
hold off;
xlim([fft_center - zoom_radius_pixels, fft_center + zoom_radius_pixels]);
ylim([fft_center - zoom_radius_pixels, fft_center + zoom_radius_pixels]);
legend('show', 'Location', 'southoutside', 'FontSize', 8);
axis on; % Show axes for clarity of zoom
xticks([]); yticks([]); % But hide tick values

% --- FFT of Illumination Pattern ---
h_ax2 = subplot(1, 3, 2);
imagesc(fft_illumination);
colormap(h_ax2, 'parula'); axis image off;
title({'FFT of Illumination Pattern';'(Zoomed View)'}, 'FontSize', 10);
hold on;
plot(fft_center - freq_illumination*cos(angle_illumination_rad), fft_center - freq_illumination*sin(angle_illumination_rad), ...
     'go', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', sprintf('Illumination Freq (%.1f)', freq_illumination));
plot(fft_center + freq_illumination*cos(angle_illumination_rad), fft_center + freq_illumination*sin(angle_illumination_rad), ...
     'go', 'MarkerSize', 15, 'LineWidth', 1.5, 'HandleVisibility','off');
hold off;
xlim([fft_center - zoom_radius_pixels, fft_center + zoom_radius_pixels]);
ylim([fft_center - zoom_radius_pixels, fft_center + zoom_radius_pixels]);
legend('show', 'Location', 'southoutside', 'FontSize', 8);
axis on; 
xticks([]); yticks([]);

% --- FFT of Moiré Image ---
h_ax3 = subplot(1, 3, 3);
imagesc(fft_moire);
colormap(h_ax3, 'parula'); axis image off;
title({'FFT of Moiré Image';'(Zoomed View)'}, 'FontSize', 10);

% Calculate frequency vector components in cycles per grid (for plotting markers)
fx_s = freq_sample * cos(angle_sample_rad);         % x-component of sample frequency vector
fy_s = freq_sample * sin(angle_sample_rad);         % y-component of sample frequency vector
fx_i = freq_illumination * cos(angle_illumination_rad); % x-component of illumination frequency vector
fy_i = freq_illumination * sin(angle_illumination_rad); % y-component of illumination frequency vector

% Moiré Difference frequency components (vector subtraction)
fx_diff = fx_s - fx_i;
fy_diff = fy_s - fy_i;
freq_diff = sqrt(fx_diff^2 + fy_diff^2); % Magnitude of the difference frequency

% Moiré Sum frequency components (vector addition)
fx_sum = fx_s + fx_i;
fy_sum = fy_s + fy_i;
freq_sum = sqrt(fx_sum^2 + fy_sum^2);   % Magnitude of the sum frequency

hold on;
% Plot original sample frequency peaks (for context within the Moiré FFT)
plot(fft_center + fx_s, fft_center - fy_s, 'ro', 'MarkerSize', 9, ...
    'DisplayName', sprintf('Sample (Orig, %.1f)', freq_sample));
plot(fft_center - fx_s, fft_center + fy_s, 'ro', 'MarkerSize', 9, 'HandleVisibility','off');

% Plot original illumination frequency peaks (for context)
plot(fft_center - fx_i, fft_center - fy_i, 'go', 'MarkerSize', 9, ...
    'DisplayName', sprintf('Illum (Orig, %.1f)', freq_illumination));
plot(fft_center + fx_i, fft_center + fy_i, 'go', 'MarkerSize', 9, 'HandleVisibility','off');

% Plot Moiré difference frequency peaks (the low-frequency component)
plot(fft_center - fx_diff, fft_center - fy_diff, 'bo',  'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Moiré Diff (Low Freq, %.1f)', freq_diff));
plot(fft_center + fx_diff, fft_center + fy_diff, 'bo',  'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility','off');

% Plot Moiré sum frequency peaks
plot(fft_center - fx_sum, fft_center - fy_sum, 'mo', 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Moiré Sum (High Freq, %.1f)', freq_sum));
plot(fft_center + fx_sum, fft_center + fy_sum, 'mo', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility','off');
hold off;

xlim([fft_center - zoom_radius_pixels*1.1, fft_center + zoom_radius_pixels*1.1]);
ylim([fft_center - zoom_radius_pixels*1.1, fft_center + zoom_radius_pixels*1.1]);
legend('show', 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);
axis on;
xticks([]); yticks([]);

sgtitle(sprintf('Zoomed FFTs: Sample Freq=%.1f, Illum Freq=%.1f. Moiré Diff Freq=%.1f', freq_sample, freq_illumination, freq_diff), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Done.\n');