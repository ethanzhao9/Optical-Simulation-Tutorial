%% 2D STED PSF simulation
% 1) Excitation PSF
% 2) STED donut from a vortex phase mask in the pupil
% 3) Effective PSF = Excitation PSF × survival probability

clear; close all; clc;
%% --- 1) Parameters ---
lambda_exc  = 488e-9;          % excitation wavelength [m]
lambda_sted = 592e-9;          % STED wavelength [m]
NA          = 1.4;             % numerical aperture
N           = 512;             % grid size
L           = 2.0e-6;          % field of view in image plane [m]

I_sat   = 1.0;                 % saturation intensity (arb.)
S0_list = [0 5 20];            % S0 = I_peak / I_sat

% Survival model 
survival_mode = 'pulsed';

%% --- 2) Excitation PSF (normalized) ---
[x, I_exc] = generate_psf(lambda_exc, NA, N, L);
I_exc = I_exc / max(I_exc(:));

%% --- 3) STED donut PSF from vortex phase mask ---
charge = 1;
I_donut = generate_vortex_donut(lambda_sted, NA, N, L, charge);
I_donut = I_donut / max(I_donut(:));

%% --- 4) Effective PSF for different STED strengths ---
num_cases   = numel(S0_list);
center_idx  = floor(N/2) + 1;
profiles    = zeros(num_cases, N);
fwhm_values = zeros(num_cases, 1);
legends     = cell(num_cases, 1);

for k = 1:num_cases
    S0     = S0_list(k);
    I_sted = (S0 * I_sat) * I_donut;

    eta   = survival_probability(I_sted, I_sat, survival_mode);
    I_eff = I_exc .* eta;
    I_eff = I_eff / max(I_eff(:));

    profiles(k,:)  = I_eff(center_idx,:);
    fwhm_values(k) = fwhm_simple(x, profiles(k,:));
    legends{k}     = sprintf('S_0=%g (FWHM=%.0f nm)', S0, fwhm_values(k)*1e9);

    if k == num_cases
        I_eff_max = I_eff;          % store strongest STED effective PSF
        eta_max   = eta;            % store strongest STED survival map
    end
end

figure('Color','w','Position',[160 90 1050 720]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% Top row: excitation / survival / effective (strongest)
nexttile(1); imagesc(x*1e9, x*1e9, I_exc);     axis image; colormap(gca,'hot');  colorbar;
title('Excitation PSF'); xlabel('x (nm)'); ylabel('y (nm)'); xlim([-400 400]); ylim([-400 400]);

nexttile(2); imagesc(x*1e9, x*1e9, eta_max);   axis image; colormap(gca,'gray'); colorbar;
title(sprintf('Survival \\eta (S_0=%g)', S0_list(end)));
xlabel('x (nm)'); ylabel('y (nm)'); xlim([-400 400]); ylim([-400 400]);

nexttile(3); imagesc(x*1e9, x*1e9, I_eff_max); axis image; colormap(gca,'hot');  colorbar;
title(sprintf('Effective PSF (S_0=%g)', S0_list(end)));
xlabel('x (nm)'); ylabel('y (nm)'); xlim([-400 400]); ylim([-400 400]);

% Bottom row: cross-sections
ax = nexttile(4,[1 3]); hold(ax,'on'); grid(ax,'on');
colors = lines(num_cases);
for k = 1:num_cases
    plot(ax, x*1e9, profiles(k,:), 'LineWidth', 2, 'Color', colors(k,:));
end
yline(ax, 0.5, 'k:', 'Half Max');
xlabel(ax, 'Position (nm)'); ylabel(ax, 'Normalized Intensity');
title(ax, 'Effective PSF cross-sections'); legend(ax, legends, 'Location','northeast');
xlim(ax, [-300 300]); hold(ax,'off');


function [x, I] = generate_psf(lambda, NA, N, L)
% Scalar PSF intensity from a circular pupil
    dx = L / N;
    x  = (-N/2 : N/2-1) * dx;

    df = 1 / L;
    f  = (-N/2 : N/2-1) * df;
    [Fx, Fy] = meshgrid(f, f);
    Rho_f = hypot(Fx, Fy);

    f_cutoff = NA / lambda;
    Pupil = (Rho_f <= f_cutoff);

    Field = fftshift(fft2(ifftshift(Pupil)));
    I = abs(Field).^2;
end

function I = generate_vortex_donut(lambda, NA, N, L, charge)
% Donut PSF from a vortex phase mask
    df = 1 / L;
    f  = (-N/2 : N/2-1) * df;
    [Fx, Fy] = meshgrid(f, f);
    [Phi, Rho_f] = cart2pol(Fx, Fy);

    f_cutoff = NA / lambda;
    Aperture = (Rho_f <= f_cutoff);

    VortexPhase = exp(1i * charge * Phi);
    PupilField  = Aperture .* VortexPhase;

    Field = fftshift(fft2(ifftshift(PupilField)));
    I = abs(Field).^2;
end

function eta = survival_probability(I_sted, I_sat, mode)
    switch lower(mode)
        case 'pulsed'
            eta = exp(-I_sted / I_sat);
        case 'cw'
            eta = 1 ./ (1 + I_sted / I_sat);
    end
end

function w = fwhm_simple(x, profile)
    half_max = max(profile) / 2;
    idx = find(profile >= half_max);
    if isempty(idx), w = 0; return; end
    w = x(idx(end)) - x(idx(1));
end
