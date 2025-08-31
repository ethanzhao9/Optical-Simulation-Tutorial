% Simulate imaging results from a 3D object with and without a pinhole
% using effective PSFs in confocal microscopy.
% Author: Ethan Zhao

clear; clc; close all;

%% ====== Optical and Simulation Parameters ======
lambda_exc = 488e-9;          % Excitation wavelength [m]
lambda_em  = 520e-9;          % Emission wavelength [m]
NA         = 0.8;             % Numerical aperture
n_medium   = 1.33;            % Refractive index of immersion medium

grid_size = 128;              % Image grid (pixels)
dx = 0.05e-6;                 % Pixel size [m]
z_range = [-5, 5]*1e-6;       % Z-scan range [m]
dz = 0.1e-6;                  % Z-step [m]
z_vals = z_range(1):dz:z_range(2);
Nz = numel(z_vals);

%% ====== Compute Excitation and Detection PSFs ======
fprintf('Computing excitation PSF...\n');
psf_exc = compute_psf(lambda_exc, NA, n_medium, grid_size, dx, z_vals);

fprintf('Computing detection PSF (no pinhole)...\n');
psf_det_raw = compute_psf(lambda_em, NA, n_medium, grid_size, dx, z_vals);

%% ====== Apply Pinhole to Detection PSF ======
airy_radius = 0.61 * lambda_em / NA;          % Airy disk radius [m]
pinhole_radius_m = 1 * airy_radius;           % Set pinhole diameter to 1 AU
pinhole_radius_px = pinhole_radius_m / dx;

x = (-grid_size/2 : grid_size/2 - 1) * dx;
[X, Y] = meshgrid(x, x);
R = sqrt(X.^2 + Y.^2);
pinhole_mask = double(R <= pinhole_radius_m); % Binary circular mask

psf_det_pinhole = zeros(size(psf_det_raw));
for zi = 1:Nz
    psf_det_pinhole(:,:,zi) = conv2(psf_det_raw(:,:,zi), pinhole_mask, 'same');
end

%% ====== Calculate Effective PSFs ======
psf_eff_raw     = psf_exc;
psf_eff_pinhole = psf_exc .* psf_det_pinhole;

% Normalize both for fair comparison
psf_eff_raw     = psf_eff_raw / max(psf_eff_raw(:));
psf_eff_pinhole = psf_eff_pinhole / max(psf_eff_pinhole(:));

%% ====== Generate a Random 3D Object ======
rng(0);                      % For reproducibility
num_points = 200;            % Number of bright spots
object = zeros(grid_size, grid_size, Nz);

for i = 1:num_points
    xi = randi(grid_size);
    yi = randi(grid_size);
    zi = randi(Nz);
    object(xi, yi, zi) = 1;
end

%% ====== Image Formation and Noise Simulation ======
fprintf('Forming and adding noise to simulated images...\n');

% 1. Convolve object with effective PSFs
img_sim_raw     = convn(object, psf_eff_raw, 'same');
img_sim_pinhole = convn(object, psf_eff_pinhole, 'same');

% 2. Scale to photon counts (simulate laser intensity and fluorophore brightness)
max_photon_count = 500;  % Peak signal level (typical for confocal)
img_raw_scaled     = img_sim_raw * max_photon_count;
img_pinhole_scaled = img_sim_pinhole * max_photon_count;

% 3. Add Poisson (photon shot) noise
img_raw_poisson     = poissrnd(img_raw_scaled);
img_pinhole_poisson = poissrnd(img_pinhole_scaled);

% 4. Add Gaussian readout noise (camera/PMT noise)
read_noise_std = 8;  % Typical value for PMT/sCMOS in photons
img_no_pinhole   = double(img_raw_poisson) + read_noise_std * randn(size(img_raw_poisson));
img_with_pinhole = double(img_pinhole_poisson) + read_noise_std * randn(size(img_pinhole_poisson));

%% ====== Extract Central Slices for Display ======
z0_idx = round(Nz / 2);
yc = round(grid_size / 2);

img_xy_raw      = img_no_pinhole(:, :, z0_idx);
img_xy_pinhole  = img_with_pinhole(:, :, z0_idx);

img_xz_raw      = squeeze(img_no_pinhole(:, yc, :))';
img_xz_pinhole  = squeeze(img_with_pinhole(:, yc, :))';

%% ====== Visualization ======
figure('Position', [100, 100, 1000, 400]); colormap hot;

subplot(2, 2, 1); imagesc(img_xy_raw); axis image off;
title('XY Slice @ z = 0 (No Pinhole)');

subplot(2, 2, 2); imagesc(img_xy_pinhole); axis image off;
title('XY Slice @ z = 0 (With Pinhole)');

subplot(2, 2, 3); imagesc(img_xz_raw); axis image off;
title('XZ Profile (No Pinhole)');

subplot(2, 2, 4); imagesc(img_xz_pinhole); axis image off;
title('XZ Profile (With Pinhole)');

%% ====== PSF Computation Function ======
function psf = compute_psf(lambda, NA, n, N, dx, z_vals)
    k = 2 * pi * n / lambda;                           % Wave number
    fx = (-N/2 : N/2 - 1) / (N * dx);                  % Spatial frequency
    [FX, FY] = meshgrid(fx, fx);
    kx = 2 * pi * FX;
    ky = 2 * pi * FY;
    k_rho = sqrt(kx.^2 + ky.^2);
    kz = real(sqrt(k^2 - k_rho.^2));                   % Propagation term

    pupil = k_rho <= (2 * pi * NA / lambda);           % Circular aperture

    Nz = numel(z_vals);
    E = zeros(N, N, Nz);

    for zi = 1:Nz
        z = z_vals(zi);
        phase = exp(1i * kz * z);
        E(:,:,zi) = fftshift(ifft2(ifftshift(pupil .* phase)));
    end

    psf = abs(E).^2;                                    % Intensity PSF
end