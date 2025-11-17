%% Filtered Backprojection (FBP) demo
% Compare UBP, manual Ram-Lak FBP, and Shepp-Logan FBP
clear; close all; clc;

%% Parameters
N     = 512;               % reconstruction image size
theta = 0:1:179;           % projection angles (degrees)
P     = phantom('Modified Shepp-Logan', N);

%% Radon transform (sinogram)
[R, xp] = radon(P, theta);           % R: [numDetectors x numAngles]
[numDet, K] = size(R);

%% 1) Unfiltered backprojection (UBP) via MATLAB iradon
UBP_iradon = iradon(R, theta, 'linear', 'none', 1, N);

%% 2) Manual FBP with Ram-Lak filter
% (a) Build 1D Ram-Lak filter in frequency domain (FFT convention)
L = 2^nextpow2(2*numDet);           % zero-padding length (>= 2*numDet)

% Frequency axis in fftshifted order: approximately [-1, 1)
k            = (-L/2 : L/2-1)' / (L/2);
ramlak_shift = abs(k);              % |ω| in shifted frequency domain
ramlak       = ifftshift(ramlak_shift);   % align with FFT output order

% (b) Filter each projection in frequency domain
R_filt  = zeros(numDet, K);         % filtered projections
midL    = floor(L/2) + 1;           % center index in padded vector
halfDet = floor(numDet/2);

for kAng = 1:K
    p = R(:, kAng);                 % projection at angle theta(kAng)

    % Zero-pad projection to length L and center it in p_pad
    p_pad = zeros(L, 1);
    idx_start = midL - halfDet;
    idx_end   = idx_start + numDet - 1;
    p_pad(idx_start:idx_end) = p;

    % FFT -> multiply by Ram-Lak filter -> IFFT
    P_fft      = fft(p_pad);
    P_fft_filt = P_fft .* ramlak;
    p_filt_pad = real(ifft(P_fft_filt));

    % Take central numDet samples back
    R_filt(:, kAng) = p_filt_pad(idx_start:idx_end);
end

% (c) Backprojection of the filtered projections (manual FBP)
[xGrid, yGrid] = meshgrid(((1:N)-(N+1)/2), -((1:N)-(N+1)/2));  % image grid (y axis up)
FBP_RamLak = zeros(N);
dtheta     = pi / K;                % numerical integration step Δθ

for kAng = 1:K
    th = deg2rad(theta(kAng));

    % Detector coordinate for each pixel: s = x cosθ + y sinθ
    s = xGrid.*cos(th) + yGrid.*sin(th);

    % Interpolate filtered projection at positions s
    proj_filt = R_filt(:, kAng);
    slice_k   = interp1(xp, proj_filt, s, 'linear', 0);

    % Accumulate over angles (discrete approximation of ∫ dθ)
    FBP_RamLak = FBP_RamLak + slice_k * dtheta;
end

%% 3) Shepp-Logan FBP via MATLAB iradon
FBP_SheppLogan = iradon(R, theta, 'linear', 'Shepp-Logan', 1, N);

%% Visualization
figure('Color','w','Position',[50 50 1600 400]); colormap gray;

subplot(1,4,1); imshow(P, []); 
title('Ground Truth (Phantom)');

subplot(1,4,2); imshow(UBP_iradon, []); 
title('UBP (iradon, ''none'')');

subplot(1,4,3); imshow(FBP_RamLak, []); 
title('Manual FBP (Ram-Lak)');

subplot(1,4,4); imshow(FBP_SheppLogan, []); 
title('FBP (iradon, Shepp-Logan)');
