% Simulates the 3D Point Spread Function (PSF) using exact angular spectrum propagation.
% Models the propagation of light from the pupil plane to multiple axial planes.
% Author: Ethan Zhao

clear; close all; clc;

%% ===================== Parameters =====================
lambda = 488e-9;           % Wavelength [m]
NA = 0.8;                  % Numerical Aperture
n_medium = 1.33;           % Refractive index of the immersion medium (e.g., water)
k = 2 * pi * n_medium / lambda;  % Wavenumber [1/m]

grid_size = 256;           % Number of pixels along one axis
dx = 0.05e-6;              % Sampling interval in the lateral (x, y) dimensions [m/pixel]

z_range = [-10, 10] * 1e-6;  % Axial simulation range [m]
dz = 0.05e-6;               % Sampling interval along z [m]
z_vals = z_range(1):dz:z_range(2);  
Nz = numel(z_vals);         % Number of z-slices

%% ===================== Frequency Coordinates =====================
% Create spatial frequency grid (Fourier domain of the image plane)
fx = (-grid_size/2 : grid_size/2 - 1) / (grid_size * dx);  
[fx_grid, fy_grid] = meshgrid(fx, fx);
kx = 2 * pi * fx_grid;       % Angular spatial frequency [rad/m]
ky = 2 * pi * fy_grid;
k_rho = sqrt(kx.^2 + ky.^2); % Radial frequency [rad/m]

% Compute z-component of wavevector, filter out evanescent components
kz = sqrt(k^2 - k_rho.^2);
kz = real(kz);  % Only include propagating waves

%% ===================== Pupil Function =====================
f_cutoff = NA / lambda;                          % Cutoff spatial frequency [1/m]
pupil_function = k_rho <= (2 * pi * f_cutoff);   % Circular aperture (ideal binary mask)

%% ===================== Field Propagation Loop =====================
E = zeros(grid_size, grid_size, Nz);  % 3D complex field volume

for zi = 1:Nz
    z = z_vals(zi);

    % Apply propagation phase to each frequency component
    defocus_phase = exp(1i * kz * z);

    % Multiply pupil by phase term
    pupil_z = pupil_function .* defocus_phase;

    % Inverse Fourier transform to compute field at axial depth z
    field_z = fftshift(ifft2(ifftshift(pupil_z)));

    % Store complex field
    E(:, :, zi) = field_z;
end

%% ===================== PSF Volume Calculation =====================
psf_volume = abs(E).^2;                   % PSF = |Field|^2
psf_volume = psf_volume / max(psf_volume(:));  % Normalize for visualization

%% ===================== Visualization =====================
% Extract center slices for inspection
[~, focus_idx] = min(abs(z_vals));        % Index of focal plane (z = 0)
y_center = round(grid_size / 2);          % Central Y slice

psf_xy = psf_volume(:, :, focus_idx);     % XY slice at focal plane
psf_xz = squeeze(psf_volume(:, y_center, :))';  % XZ slice through center Y

% Plot
figure('Position', [100, 100, 1000, 420]);

subplot(1, 2, 1);
imagesc(psf_xy);
axis image off;
colormap hot;
title('Simulated 3D PSF (XY Focal Plane)');

subplot(1, 2, 2);
imagesc(psf_xz); 
axis image off;
title('XZ Plane (Axial Cross Section)');
colormap hot;
title('Simulated 3D PSF (XZ Cross Section)');

ax = gcf;
exportgraphics(ax,'PSF_3D.jpg','Resolution',300);