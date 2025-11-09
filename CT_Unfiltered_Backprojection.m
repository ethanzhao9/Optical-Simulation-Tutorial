%% Unfiltered Backprojection (UBP) demo
clear; close all; clc;

%% Parameters
N     = 512;            % image size
theta = 0:1:179;        % projection angles (denser -> better UBP)

P = phantom('Modified Shepp-Logan', N);

%% Radon transform (sinogram)
[R, xp] = radon(P, theta);              % R: [numDetectors x numAngles]
K       = numel(theta);

%% Manual Unfiltered Backprojection (UBP)
% Build a Cartesian grid centered at the image origin.
[xGrid, yGrid] = meshgrid(((1:N)-(N+1)/2), -((1:N)-(N+1)/2));

B = zeros(N);                           % accumulator
dtheta = pi / K;                        % numerical integration step

for k = 1:K
    th = deg2rad(theta(k));

    % For each pixel, compute detector coordinate: s = x cosθ + y sinθ
    s = xGrid.*cos(th) + yGrid.*sin(th);

    % Interpolate the 1D projection R(:,k) at positions s and "smear back"
    slice_k = interp1(xp, R(:,k), s, 'linear', 0);

    % Accumulate (discrete approximation to the θ integral)
    B = B + slice_k * dtheta;
end

UBP_manual = B;

%% MATLAB built-in unfiltered backprojection (for comparison)
UBP_iradon = iradon(R, theta, 'linear', 'none', 1, N);  % no filter => UBP

%% Visualization
figure('Color','w','Position',[80 80 1200 400]); colormap gray;

subplot(1,3,1); imshow(P,[]);        title('Ground Truth (Phantom)');
subplot(1,3,2); imshow(UBP_manual,[]); title('Manual Unfiltered Backprojection (UBP)');
subplot(1,3,3); imshow(UBP_iradon,[]); title('MATLAB iradon (''none'') = UBP');
