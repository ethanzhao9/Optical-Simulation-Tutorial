%% Demo for phase retrieval: Creating 3x3-Spot Hologram 
%% Author: Ethan Zhao
clear; clc; close all;

%% 1. Simulation Parameters
lambda      = 488;        % Wavelength (nm)
NA          = 1.0;         % Numerical aperture
n_medium    = 1.0;         % Refractive index
N           = 256;         % Grid size (pixels)
dxy         = 150;         % Pupil plane pixel size (nm)
nIter       = 150;         % Number of GS iterations
alpha       = 0.9;         % Damping factor

%% 2. Pupil and Coordinates
R_pupil = NA * N * dxy / (n_medium * lambda);
[x, y] = meshgrid(-(N/2):(N/2-1));
r = sqrt(x.^2 + y.^2);
pupilAmp = double(r <= R_pupil);

%% 3. Target Amplitude (3x3 Gaussian Spots)
gridSize = 3;
sep = 50;
targetAmp = zeros(N);
[xg, yg] = meshgrid(-(gridSize-1)/2:(gridSize-1)/2);
spotPos = [xg(:), yg(:)] * sep;

for k = 1:size(spotPos,1)
    cx = round(spotPos(k,1) + N/2 + 1);
    cy = round(spotPos(k,2) + N/2 + 1);
    targetAmp(cy, cx) = 1;
end

% Apply Gaussian blur
h = fspecial('gaussian', 7, 1.5);
targetAmp = imfilter(targetAmp, h, 'replicate');

% Normalize
targetAmp = targetAmp / sqrt(sum(targetAmp(:).^2));
signalMask = targetAmp > 0.1 * max(targetAmp(:));

%% 4. Gerchberg-Saxton Algorithm
pupilPhase = 2 * pi * rand(N);
pupilField = pupilAmp .* exp(1j * pupilPhase);
effHist = zeros(1, nIter);
uniHist = zeros(1, nIter);

for iter = 1:nIter
    fieldFocal = fftshift(fft2(ifftshift(pupilField)));
    I_focal = abs(fieldFocal).^2;

    % Metrics
    effHist(iter) = sum(I_focal(signalMask)) / sum(I_focal(:));
    peaks = zeros(1, size(spotPos,1)); 
    w = 3;
    for k = 1:size(spotPos,1)
        cx = round(spotPos(k,1) + N/2 + 1);
        cy = round(spotPos(k,2) + N/2 + 1);
        region = I_focal(cy-w:cy+w, cx-w:cx+w);
        peaks(k) = max(region(:));
    end
    uniHist(iter) = 1 - (max(peaks) - min(peaks)) / (max(peaks) + min(peaks) + eps);

    % Constraint update
    fieldFocal = targetAmp .* exp(1j * angle(fieldFocal));
    fieldBack = ifftshift(ifft2(fftshift(fieldFocal)));
    pupilField = (1 - alpha) * pupilField + alpha * pupilAmp .* exp(1j * angle(fieldBack));
end

%% 5. Results
I_final = I_focal / max(I_focal(:));
figure; imagesc(I_final); axis image; colormap hot; colorbar;
title('Focal Plane Intensity'); xlabel('x'); ylabel('y');

figure; imagesc(angle(pupilField) .* pupilAmp); axis image;
colormap hsv; colorbar; title('Pupil Phase Mask');

figure;
subplot(1,2,1); plot(1:nIter, effHist*100, 'LineWidth', 1.5);
title('Efficiency'); xlabel('Iteration'); ylabel('%'); grid on; ylim([0 100]);
subplot(1,2,2); plot(1:nIter, uniHist*100, 'LineWidth', 1.5);
title('Uniformity'); xlabel('Iteration'); ylabel('%'); grid on; ylim([0 100]);

fprintf('\nFinal Efficiency: %.2f%%\n', effHist(end)*100);
fprintf('Final Uniformity: %.2f%%\n', uniHist(end)*100);
