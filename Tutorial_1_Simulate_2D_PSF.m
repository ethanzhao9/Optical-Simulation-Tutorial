% Simulates the 2D Point Spread Function (PSF) using Fourier optics principles.
% Models the transformation from the pupil plane (angular spectrum) to the image (focal) plane.
% Author: Ethan Zhao

clear; close all; clc;

%% ===================== Parameters =====================
lambda = 488e-9;        % Wavelength [m] (488 nm, blue light)
NA = 0.2;               % Numerical Aperture of the objective lens
n_medium = 1.33;        % Refractive index of the immersion medium (e.g., water)
k = 2 * pi / lambda;    % Wavenumber [1/m]

grid_size = 512;        % Number of pixels along one axis
dx_image = 0.05e-6;     % Sampling interval in the image (focal) plane [m/pixel]

% Spatial frequency coordinates (Fourier domain of image plane)
fx = (-grid_size/2 : grid_size/2 - 1) / (grid_size * dx_image);  
[fx_grid, fy_grid] = meshgrid(fx, fx);
rho = sqrt(fx_grid.^2 + fy_grid.^2);   % Radial spatial frequency [1/m]

%% ===================== Pupil Function =====================
% The pupil function defines the angular acceptance of the imaging system.
% A circular aperture with cutoff frequency determined by NA
f_cutoff = NA / lambda;                     % Maximum spatial frequency (diffraction limit)
pupil_function = rho <= f_cutoff;          % Circular binary mask (ideal pupil)

%% ===================== Field & PSF Calculation =====================
% The focal plane field is the Fourier transform of the pupil function
focal_plane_field = ifftshift(ifft2(fftshift(pupil_function)));

% PSF is the intensity (magnitude squared) of the complex field
psf_intensity = abs(focal_plane_field).^2;

% Normalize the intensity for visualization
psf_intensity = psf_intensity / max(psf_intensity(:));

%% ===================== Visualization =====================
figure;
imagesc(psf_intensity);
axis image off;
colorbar;
title('Simulated 2D PSF at the Focal Plane');



