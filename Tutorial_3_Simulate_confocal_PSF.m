% Simulates the 3D Point Spread Function (PSF) of a confocal microscope using angular spectrum propagation.
% Author: Ethan Zhao 
clear; clc; close all;

%% ===================== Optical Parameters =====================
lambda_exc = 488e-9;        % Excitation wavelength [m]
lambda_em  = 520e-9;        % Emission wavelength [m]
NA         = 0.8;           % Numerical Aperture
n_medium   = 1.33;          % Refractive index of medium

%% ===================== Spatial Sampling Parameters =====================
grid_size = 256;              % Grid size (pixels)
dx = 0.05e-6;                 % Lateral sampling [m/pixel]
z_range = [-5, 5]*1e-6;       % Axial range [m]
dz = 0.05e-6;                 % Axial step [m]
z_vals = z_range(1):dz:z_range(2);
Nz = numel(z_vals);

%% ===================== Pinhole Definition =====================
pinhole_AU = 1.0;                                 % Diameter in Airy Units
airy_radius = 0.61 * lambda_em / NA;              % Airy disk radius [m]
pinhole_radius = (pinhole_AU / 2) * airy_radius;  % Pinhole radius [m]
pinhole_radius_px = pinhole_radius / dx;

%% ===================== Excitation PSF =====================
psf_exc = compute_PSF(lambda_exc, NA, n_medium, grid_size, dx, z_vals);

%% ===================== Detection PSF =====================
psf_det = compute_PSF(lambda_em, NA, n_medium, grid_size, dx, z_vals);

%% ===================== Apply Pinhole to Detection PSF =====================
x = (-grid_size/2 : grid_size/2 - 1) * dx;
[X, Y] = meshgrid(x, x);
R = sqrt(X.^2 + Y.^2);
pinhole_mask = double(R <= pinhole_radius);
psf_det_filtered = zeros(size(psf_det));

for zi = 1:Nz
    psf_det_filtered(:,:,zi) = conv2(psf_det(:,:,zi), pinhole_mask, 'same');
end

%% ===================== Confocal PSF =====================
psf_confocal = psf_exc .* psf_det_filtered;
psf_confocal = psf_confocal / max(psf_confocal(:));

psf_noPinhole = psf_exc .* psf_det;
psf_noPinhole = psf_noPinhole / max(psf_noPinhole(:));
%% ===================== Extract Slices =====================
[~, z0_idx] = min(abs(z_vals));
y_center = round(grid_size / 2);

psf_exc_xy = psf_exc(:,:,z0_idx);
psf_exc_xz = squeeze(psf_exc(:, y_center, :))';

psf_conf_xy = psf_confocal(:,:,z0_idx);
psf_conf_xz = squeeze(psf_confocal(:, y_center, :))';

psf_noPinhole_xy = psf_noPinhole(:,:,z0_idx);
psf_noPinhole_xz = squeeze(psf_noPinhole(:, y_center, :))';

psf_exc_xy = psf_exc_xy / max(psf_exc_xy(:));
psf_exc_xz = psf_exc_xz / max(psf_exc_xz(:));

%% ===================== Visualization =====================
figure('Position', [100, 100, 1200, 800]); colormap hot;

subplot(2, 2, 1);
imagesc(psf_exc_xy); axis image off; colorbar;
title('Excitation PSF (XY, z=0)');

subplot(2, 2, 2);
imagesc(psf_exc_xz); axis image off; colorbar;
title('Excitation PSF (XZ)');

subplot(2, 2, 3);
imagesc(psf_conf_xy); axis image off; colorbar;
title('Confocal PSF (XY)');

subplot(2, 2, 4);
imagesc(psf_conf_xz); axis image off; colorbar;
title('Confocal PSF (XZ)');

% ax = gcf;
% exportgraphics(ax,'Confocal_PSF.jpg','Resolution',300);
%% ===================== Profile Comparison =====================
x_um = x * 1e6;
z_um = z_vals * 1e6;

figure('Position', [43         329        1723         400]);

subplot(1, 2, 1);
plot(x_um, psf_exc_xy(y_center,:), 'b', 'LineWidth', 1.5); hold on;
plot(x_um, psf_conf_xy(y_center,:), 'r--', 'LineWidth', 1.5);
plot(x_um, psf_noPinhole_xy(y_center,:), 'k--', 'LineWidth', 1.5);
xlabel('X (\mum)'); ylabel('Intensity');
title('Lateral Profile (z=0)');
legend('Excitation', 'Confocal', 'No Pinhole'); grid on;
xlim([-2,2]);

subplot(1, 2, 2);
plot(z_um, psf_exc_xz(:, y_center), 'b', 'LineWidth', 1.5); hold on;
plot(z_um, psf_conf_xz(:, y_center), 'r--', 'LineWidth', 1.5);
plot(z_um, psf_noPinhole_xz(:, y_center), 'k--', 'LineWidth', 1.5);
xlabel('Z (\mum)'); ylabel('Intensity');
title('Axial Profile');
legend('Excitation', 'Confocal','No Pinhole'); grid on;
xlim([-5,5]);

% ax = gcf;
% exportgraphics(ax,'Confocal_PSF_profile.jpg','Resolution',300);

function psf = compute_PSF(lambda, NA, n, N, dx, z_vals)
    k = 2 * pi * n / lambda;
    fx = (-N/2 : N/2 - 1) / (N * dx);
    [FX, FY] = meshgrid(fx, fx);
    kx = 2*pi*FX; ky = 2*pi*FY;
    k_rho = sqrt(kx.^2 + ky.^2);
    kz = real(sqrt(k^2 - k_rho.^2));
    pupil = k_rho <= (2*pi*NA/lambda);

    Nz = numel(z_vals);
    E = zeros(N, N, Nz);
    for zi = 1:Nz
        phase = exp(1i * kz * z_vals(zi));
        E(:,:,zi) = fftshift(ifft2(ifftshift(pupil .* phase)));
    end
    psf = abs(E).^2;
end
