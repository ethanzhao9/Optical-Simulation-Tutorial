% Simulates the effect of Zernike aberrations on the 2D Point Spread Function (PSF).
% Models a wavefront distortion in the pupil plane and its effect on the focal plane image.
% Author: Ethan Zhao

clear; close all; clc;

%% ===================== Parameters =====================
lambda = 488e-9;            % Wavelength [m] (e.g., green light)
NA = 0.8;                   % Numerical Aperture of the objective lens
n_medium = 1.0;             % Refractive index of the imaging medium (e.g., air)

grid_size = 512;            % Number of pixels along one axis
dx_image = 50e-9;           % Sampling interval in the image plane [m/pixel]

% Spatial frequency coordinates, corresponding to the Fourier domain of the image plane
fx = (-grid_size/2 : grid_size/2 - 1) / (grid_size * dx_image);
[fx_grid, fy_grid] = meshgrid(fx, fx);
rho = sqrt(fx_grid.^2 + fy_grid.^2); % Radial spatial frequency [1/m]

%% ===================== Pupil Coordinates & Function =====================
% Define a normalized coordinate system (rho_pupil, theta_pupil) over the pupil,
% where rho_pupil = 1 at the edge of the aperture.
f_cutoff = NA / lambda;             % Maximum spatial frequency transmitted by the lens
rho_pupil = rho / f_cutoff;         % Normalized radial pupil coordinate
theta_pupil = atan2(fy_grid, fx_grid); % Azimuthal pupil coordinate

% The pupil function is a binary mask defining the system's aperture
pupil_function = rho_pupil <= 1;

%% ===================== Zernike Polynomial Definitions =====================
% Zernike polynomials (ANSI standard) are defined on the normalized pupil grid.
% The aberration strength is the RMS wavefront error in radians.
aberration_rms_rad = pi; % A value of pi corresponds to an RMS error of lambda/2.

Z_tilt        = 2 * rho_pupil .* cos(theta_pupil);
Z_defocus     = sqrt(3) * (2 * rho_pupil.^2 - 1);
Z_astigmatism = sqrt(6) * rho_pupil.^2 .* cos(2 * theta_pupil);
Z_coma        = sqrt(8) * (3 * rho_pupil.^3 - 2 * rho_pupil) .* cos(theta_pupil);
Z_spherical   = sqrt(5) * (6 * rho_pupil.^4 - 6 * rho_pupil.^2 + 1);

% List of aberrations to simulate
AberrationList = {
    zeros(grid_size, grid_size), 'No Aberration';
    Z_tilt,                      'Tilt (Z_1^1)';
    Z_defocus,                   'Defocus (Z_2^0)';
    Z_astigmatism,               'Astigmatism (Z_2^2)';
    Z_coma,                      'Coma (Z_3^1)';
    Z_spherical,                 'Spherical (Z_4^0)';
};

%% ===================== PSF Simulation & Visualization =====================
num_aberrations = size(AberrationList, 1);
figure('Position', [100, 100, 1400, 550], 'Color', 'w'); % Adjusted for new layout

% Create a tiled layout with compact spacing
t = tiledlayout(2, num_aberrations, 'TileSpacing', 'compact', 'Padding', 'compact');

ROI = 80; % Half-width for zoomed PSF display window [pixels]
cx = grid_size/2 + 1; cy = grid_size/2 + 1;

for i = 1:num_aberrations
    aberration_profile = AberrationList{i, 1};
    label = AberrationList{i, 2};

    % The complex pupil field includes the phase aberration
    wavefront_phase = aberration_rms_rad * aberration_profile;
    pupil_field = pupil_function .* exp(1j * wavefront_phase);

    % The focal plane field is the Fourier transform of the pupil field
    focal_plane_field = ifftshift(ifft2(fftshift(pupil_field)));

    psf_intensity = abs(focal_plane_field).^2;
    psf_intensity = psf_intensity / max(psf_intensity(:));

    % --- Plot Pupil Phase 
    nexttile;
    imagesc(wavefront_phase .* pupil_function);
    axis image off;
    colormap(gca, 'hsv'); colorbar;
    title(label,'FontSize',16);

    % --- Plot Resulting PSF 
    nexttile(i + num_aberrations); % Replaces subplot
    imagesc(psf_intensity(cy-ROI:cy+ROI, cx-ROI:cx+ROI));
    axis image off;
    colormap(gca, 'hot'); colorbar;
    title(label,'FontSize',16);
end