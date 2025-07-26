% Simulate noise-free raw SIM images with realistic optical parameters.
% Author: Ethan Zhao
clear; clc; close all;

%% -- Parameters ----------------------------------------------------------
N              = 256;         % Image size (pixels)
dx             = 40e-9;       % Pixel size [m]
lambda_exc     = 488e-9;      % Excitation wavelength [m]
lambda_em      = 520e-9;      % Emission wavelength [m]
NA             = 1.20;        % Numerical aperture
n_medium       = 1.33;        % Refractive index

num_orients    = 3;                              % Number of orientations
num_phases     = 3;                              % Number of phase shifts
orient_deg     = (0:num_orients-1) * 180 / num_orients;
phase_rad      = (0:num_phases-1) * 2*pi / num_phases;
pattern_factor = 0.98;                           % Fraction of max frequency
mod_depth      = 0.90;                           % Global modulation depth

%% -- Detection PSF (2D, z=0) ---------------------------------------------
psf_det = compute_PSF(lambda_em, NA, n_medium, N, dx, 0);
psf_det = psf_det / sum(psf_det(:));             % Normalized energy

%% -- Ground Truth Object -------------------------------------------------
rng(0);
obj = zeros(N);
[Xc, Yc] = meshgrid((1:N) - N/2 - 0.5);
sigma_px = 40e-9 / dx;

for i = 1:600
    cx = (rand - 0.5) * N;
    cy = (rand - 0.5) * N;
    amp = max(0, 1 + 0.2 * randn);
    obj = obj + amp * exp(-((Xc - cx).^2 + (Yc - cy).^2) / (2 * sigma_px^2));
end
obj = obj / max(obj(:));

%% -- Illumination Patterns -----------------------------------------------
k_max   = (2 * NA) / lambda_exc;
k_illum = pattern_factor * k_max;
k_pix   = k_illum * dx;

[Xp, Yp] = meshgrid((0:N-1) - N/2);
illum = zeros(N, N, num_phases, num_orients);

for o = 1:num_orients
    angle = deg2rad(orient_deg(o));
    arg = 2 * pi * k_pix * (Xp * cos(angle) + Yp * sin(angle));
    for p = 1:num_phases
        phi = phase_rad(p);
        illum(:, :, p, o) = 0.5 * (1 + mod_depth * cos(arg + phi));
    end
end

%% -- Simulate Raw SIM Images ---------------------------------------------
raw = zeros(N, N, num_phases, num_orients);
OTF_det = fft2(ifftshift(psf_det));

for o = 1:num_orients
    for p = 1:num_phases
        img_modulated = obj .* illum(:, :, p, o);
        raw(:, :, p, o) = real(ifft2(fft2(img_modulated) .* OTF_det));
    end
end

%% -- Visualization -------------------------------------------------------
figure('Position', [60, 80, 1200, 700]); 

subplot(2, 3, 1);
imagesc(obj); axis image off; 
title('Ground Truth Object');

subplot(2, 3, 2);
imagesc(illum(:, :, 1, 1)); axis image off; colormap gray;
title(sprintf('Illum. 0° / Phase 0°  (m=%.2f)', mod_depth));

subplot(2, 3, 3);
imagesc(raw(:, :, 1, 1)); axis image off;
title('Raw SIM Image (0°, 0°)');

for p = 1:num_phases
    subplot(2, 3, 3 + p);
    imagesc(raw(:, :, p, 1)); axis image off;
    title(sprintf('Raw Phase %d°', round(rad2deg(phase_rad(p)))));
end


%% OPTIONAL:  save stack as 16‑bit TIFF
stack = reshape(raw,N,N,[]);
stack = uint16(stack/max(stack(:))*65535);
imstackwrite(stack,'SIM_raw.tif');   

stack = reshape(illum,N,N,[]); 
stack = uint16(stack/max(stack(:))*65535);
imstackwrite(stack,'SIM_raw_Illumination.tif');  

function imstackwrite(img, fname)
    imwrite(img(:,:,1),fname);
    for i = 2:size(img,3), imwrite(img(:,:,i),fname,'WriteMode','append'); end
end

% -- Local Function: PSF Generator ---------------------------------------
function psf2d = compute_PSF(lambda, NA, n, N, dx, z)
    k = 2 * pi * n / lambda;
    fx = (-N/2:N/2-1) / (N * dx);
    [FX, FY] = meshgrid(fx);
    k_rho = 2 * pi * sqrt(FX.^2 + FY.^2);
    kz = real(sqrt(k^2 - k_rho.^2));
    pupil = k_rho <= 2 * pi * NA / lambda;
    phase = exp(1i * kz * z);
    E = fftshift(ifft2(ifftshift(pupil .* phase)));
    psf2d = abs(E).^2;
end