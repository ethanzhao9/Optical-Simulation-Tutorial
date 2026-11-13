% Simulate a 3-D zero-order Bessel-like beam from a hard annular pupil.
% Compare ring radius and ring width using XY sections and true XZ projections.

clear; clc; close all;

%% ===================== Optical Parameters =====================
lambda = 920e-9;                   % Vacuum excitation wavelength [m]
n_medium = 1.33;                   % Sample refractive index
na_ring = 0.32;                    % Reference annulus center, n*sin(theta)
dna_ring = 0.045;                  % Reference full radial width in NA units
na_radii = [0.28 na_ring 0.36];     % Annulus centers for radius scan
na_widths = [0.035 dna_ring 0.055]; % Full radial widths for width scan

%% ===================== Sampling and Output =====================
x = (-100:100)*0.05e-6;            % Transverse x coordinates [m]
y = x;                             % Transverse y coordinates [m]
z = (-200:200)*0.50e-6;            % Signed distance from the focus [m]
num_f = 801;                       % Quadrature points across the open annulus
export_figures = true;
save_volume = false;               % Optional reference volume as a MAT-file
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
output_dir = fullfile(script_dir,'figures');

[X,Y] = meshgrid(x,y);
r_xy = hypot(X,Y);
dr = 0.025e-6;                     % Fine radial sampling for interpolation [m]
r = (0:ceil(max(r_xy(:))/dr)+1)*dr;
[~, ix0] = min(abs(x));
[~, iy0] = min(abs(y));
[~, iz0] = min(abs(z));
x_um = x*1e6;
y_um = y*1e6;
z_um = z*1e6;

%% ===================== Hard Annulus and Parameter Scans =====================
f_max = n_medium/lambda;                         % Propagating-wave limit [cycles/m]
f_ref = na_ring/lambda;
delta_f_ref = dna_ring/lambda;
f_centers = [na_radii, na_ring*ones(size(na_widths))]/lambda;
ring_widths = [dna_ring*ones(size(na_radii)), na_widths]/lambda;

% Keep one full 3-D volume; for other cases collect views plane by plane.
reference_idx = find(f_centers == f_ref & ring_widths == delta_f_ref,1);
if isempty(reference_idx)
    reference_idx = numel(f_centers)+1;
    f_centers(reference_idx) = f_ref;
    ring_widths(reference_idx) = delta_f_ref;
end
volume = [];
cases = cell(size(f_centers));
for ii = 1:numel(f_centers)
    duplicate = find(f_centers(1:ii-1) == f_centers(ii) & ...
        ring_widths(1:ii-1) == ring_widths(ii),1);
    if ~isempty(duplicate)
        cases{ii} = cases{duplicate};
        continue
    end
    [cases{ii}, volume_case] = simulate_annulus(r,r_xy,x,z,ix0,iy0,iz0, ...
        f_centers(ii),ring_widths(ii),f_max,num_f,ii == reference_idx);
    if ii == reference_idx, volume = volume_case; end
end
ref = cases{reference_idx};
radius_cases = cases(1:numel(na_radii));
width_cases = cases(numel(na_radii)+(1:numel(na_widths)));

fprintf('Reference: %.0f nm, n = %.2f, ring NA %.3f, full width %.3f\n', ...
    lambda*1e9,n_medium,na_ring,dna_ring);
fprintf('Radius scan: center NA | XY intensity FWHM (um) | z50 (um)\n');
for ii = 1:numel(radius_cases)
    fprintf('%6.3f %10.3f %10.2f\n',na_radii(ii), ...
        radius_cases{ii}.fwhm*1e6,radius_cases{ii}.z50*1e6);
end
fprintf('Width scan: width in NA units | XY intensity FWHM (um) | z50 (um)\n');
for ii = 1:numel(width_cases)
    fprintf('%6.3f %10.3f %10.2f\n',na_widths(ii), ...
        width_cases{ii}.fwhm*1e6,width_cases{ii}.z50*1e6);
end

%% ===================== Figure 1: 3-D Beam and Orthogonal Views =====================
fig3d = figure('Color','w','Position',[50 50 1520 870]);
layout3d = tiledlayout(fig3d,2,5,'TileSpacing','compact','Padding','compact');
title(layout3d,sprintf('Optical intensity | %.0f nm, ring NA %.3f, width %.3f', ...
    lambda*1e9,na_ring,dna_ring));

ax3d = nexttile(layout3d,1,[2 1]);
% volume already follows meshgrid convention: [y, x, z].
ix = find(abs(x_um)<=2); iy = find(abs(y_um)<=2);
surface_mesh = isosurface(x_um(ix),y_um(iy),z_um,volume(iy,ix,:),0.5);
if ~isempty(surface_mesh.vertices)
    patch(ax3d,surface_mesh,'FaceColor',[0.66 0.19 0.12],'EdgeColor','none');
end
ref_z_view = min(z_um(end),max(50,1.35*ref.z50*1e6));
if ~isfinite(ref_z_view), ref_z_view = z_um(end); end
view(ax3d,[-36 12]); grid(ax3d,'on'); box(ax3d,'on');
set(ax3d,'Projection','orthographic');
xlim(ax3d,[-4 4]); ylim(ax3d,[-4 4]); zlim(ax3d,[-ref_z_view ref_z_view]);
daspect(ax3d,[1 1 1]);
xticks(ax3d,0); yticks(ax3d,0);
xlabel(ax3d,'X (µm)'); ylabel(ax3d,'Y (µm)'); zlabel(ax3d,'Z (µm)');
title(ax3d,{'3-D half-intensity core','Equal physical scale'});
axes(ax3d); camlight headlight; lighting flat;

ax_section = nexttile(layout3d,2,[2 1]);
plot_xz(ax_section,x_um,z_um,ref.xz,ref_z_view,{'XZ center section','Y = 0'});
ax_projection = nexttile(layout3d,3,[2 1]);
plot_xz(ax_projection,x_um,z_um,ref.xz_mip,ref_z_view,{'XZ max projection','Across Y'});
if isfinite(ref.z50)
    % Small marks beside the beam show the two on-axis half-intensity crossings.
    for z_half_um = [-1 1]*ref.z50*1e6
        line(ax_projection,[1.5 3.5],[z_half_um z_half_um], ...
            'Color',[0.56 0.24 0.12],'LineWidth',1.1);
    end
end

ax = nexttile(layout3d,4);
plot_xy(ax,x_um,y_um,ref.xy,'XY at focus');
if isfinite(ref.z50)
    [~, iz_side] = min(abs(z-ref.z50));
else
    [~, iz_side] = min(abs(z-z(end)/2));
end
ax = nexttile(layout3d,5);
plot_xy(ax,x_um,y_um,volume(:,:,iz_side),sprintf('XY at Z = +%.1f µm',z_um(iz_side)));

ax_detail = nexttile(layout3d,9,[1 2]);
imagesc(ax_detail,x_um,z_um,ref.xz.',[0 1]); axis(ax_detail,'xy');
xlim(ax_detail,[-4 4]); ylim(ax_detail,[-6 6]);
daspect(ax_detail,[1 1 1]); colormap(ax_detail,hot);
xlabel(ax_detail,'X (µm)'); ylabel(ax_detail,'Z (µm)');
title(ax_detail,'XZ center section | equal-scale zoom near focus');

%% ===================== Figures 2 and 3: Pupil, XY and XZ =====================
pupil_span = 1.10*max(f_centers+ring_widths/2); % Common scale across both pupil figures
fig_radius = plot_scan(radius_cases,compose('NA %.2f',na_radii), ...
    x_um,y_um,z_um,lambda,pupil_span,sprintf('%.0f nm | vary ring radius at width %.3f',lambda*1e9,dna_ring));
fig_width = plot_scan(width_cases,compose('Width %.3f',na_widths), ...
    x_um,y_um,z_um,lambda,pupil_span,sprintf('%.0f nm | vary ring width at center NA %.2f',lambda*1e9,na_ring));

%% ===================== Export =====================
if export_figures || save_volume
    if ~isfolder(output_dir), mkdir(output_dir); end
end
if export_figures
    exportgraphics(fig3d,fullfile(output_dir,'00_beam_multiscale.png'),'Resolution',300);
    exportgraphics(fig_radius,fullfile(output_dir,'01_radius_fine_scan.png'),'Resolution',300);
    exportgraphics(fig_width,fullfile(output_dir,'02_thin_ring_scan.png'),'Resolution',300);
end
if save_volume
    % volume is ordered [y, x, z]; x, y and z are in metres.
    save(fullfile(output_dir,'bessel_volume.mat'),'volume','x','y','z');
end

%% ===================== Local Functions =====================
function [beam,volume] = simulate_annulus(r,r_xy,x,z,ix0,iy0,iz0,f0,delta_f,f_max,num_f,keep_volume)
    f_inner = f0-delta_f/2;
    f_outer = f0+delta_f/2;
    if ~(delta_f > 0 && f_inner >= 0 && f_outer < f_max)
        error('Choose a positive ring width and keep both edges inside the propagating band.');
    end

    % Binary pupil: amplitude = 1 between the edges and 0 elsewhere.
    % Integrate only over the transmitting band; there is no Gaussian weight.
    f = linspace(f_inner,f_outer,num_f);
    weights = f;                            % Polar measure f df; constant df cancels
    weights([1 end]) = weights([1 end])/2;  % Trapezoidal quadrature
    kz = 2*pi*sqrt(f_max^2-f.^2);
    kz0 = 2*pi*sqrt(f_max^2-f0^2);
    phase_z = exp(1i*(kz(:)-kz0).*z);
    radial_kernel = besselj(0,2*pi*r(:)*f);
    field = (radial_kernel.*weights)*phase_z/sum(weights);
    intensity_rz = abs(field).^2;           % Same I(0,0) reference for all z

    % Sample the 3-D intensity one XY plane at a time, in [y, x, z] order.
    volume = [];
    if keep_volume, volume = zeros([size(r_xy),numel(z)],'single'); end
    beam.xy = [];
    beam.xz = zeros(numel(x),numel(z),'single');
    beam.xz_mip = zeros(numel(x),numel(z),'single');
    radial_idx = floor(r_xy/(r(2)-r(1)))+1;
    frac = (r_xy-reshape(r(radial_idx),size(r_xy)))/(r(2)-r(1));
    for iz = 1:numel(z)
        radial_profile = intensity_rz(:,iz);
        plane = single((1-frac).*radial_profile(radial_idx) + ...
            frac.*radial_profile(radial_idx+1));
        if keep_volume, volume(:,:,iz) = plane; end
        if iz == iz0, beam.xy = plane; end
        beam.xz(:,iz) = plane(iy0,:).';            % Central slice at y = 0
        beam.xz_mip(:,iz) = max(plane,[],1).';    % Project across the y rows
    end
    beam.on_axis = double(beam.xz(ix0,:));
    beam.fwhm = 2*first_half_crossing(r,intensity_rz(:,iz0));
    beam.intensity_rz = intensity_rz;
    beam.z50 = first_half_crossing(z(iz0:end),beam.on_axis(iz0:end));
    beam.f0 = f0;
    beam.delta_f = delta_f;
end

function cross = first_half_crossing(axis_values,profile)
    j = find(profile <= 0.5,1,'first');
    cross = NaN;
    if ~isempty(j) && j > 1
        t = (0.5-profile(j-1))/(profile(j)-profile(j-1));
        cross = axis_values(j-1)+t*(axis_values(j)-axis_values(j-1));
    end
end

function fig = plot_scan(cases,labels,x_um,y_um,z_um,lambda,f_span,heading)
    n_cases = numel(cases);
    fig = figure('Color','w','Position',[40 40 1440 1040]);
    cols = n_cases+1;
    layout = tiledlayout(fig,5,cols,'TileSpacing','compact','Padding','compact');
    title(layout,heading);
    z_view = min(z_um(end),60);   % Same vertical limits for every projection
    f_axis = linspace(-f_span,f_span,801);
    [FX,FY] = meshgrid(f_axis);
    f_r = hypot(FX,FY);
    for ii = 1:n_cases
        b = cases{ii};
        ax = nexttile(layout,ii);
        pupil = abs(f_r-b.f0) <= b.delta_f/2;
        imagesc(ax,f_axis*lambda,f_axis*lambda,pupil,[0 1]);
        axis(ax,'image'); axis(ax,'xy'); colormap(ax,gray(2));
        xlabel(ax,'NA_x'); ylabel(ax,'NA_y');
        title(ax,{sprintf('Pupil | %s',labels{ii}), ...
            sprintf('Center %.3f | width %.3f',b.f0*lambda,b.delta_f*lambda)});

        ax = nexttile(layout,cols+ii);
        plot_xy(ax,x_um,y_um,b.xy,sprintf('XY | FWHM %.2f µm',b.fwhm*1e6));
        ax = nexttile(layout,2*cols+ii,[3 1]);
        plot_xz(ax,x_um,z_um,b.xz_mip,z_view, ...
            {'XZ max projection',sprintf('z_{50} = %.1f µm',b.z50*1e6)});
    end
    ax = nexttile(layout,cols,[2 1]);
    axis(ax,'off');
    text(ax,0.05,0.65,{'XZ: equal physical scale','X and Z both in µm','XY panels are zoomed'}, ...
        'FontSize',12,'Units','normalized','VerticalAlignment','top');
    ax = nexttile(layout,3*cols,[3 1]);
    profiles = cell2mat(cellfun(@(b) b.on_axis(:),cases,'UniformOutput',false));
    plot(ax,profiles,z_um, ...
        'LineWidth',1.5);
    hold(ax,'on'); xline(ax,0.5,':','HandleVisibility','off');
    xlabel(ax,'On-axis intensity'); ylabel(ax,'Z (µm)');
    xlim(ax,[0 1.05]); ylim(ax,[-z_view z_view]); grid(ax,'on');
    title(ax,'On-axis intensity vs Z');
    legend(ax,labels,'Location','northeast');
end

function plot_xy(ax,x_um,y_um,img,heading)
    imagesc(ax,x_um,y_um,img,[0 1]); axis(ax,'image'); axis(ax,'xy');
    xlim(ax,[-4 4]); ylim(ax,[-4 4]);
    colormap(ax,hot); xlabel(ax,'X (µm)'); ylabel(ax,'Y (µm)');
    title(ax,heading);
end

function plot_xz(ax,x_um,z_um,img,z_view,heading)
    % Transpose [x, z] into [z, x]: Z is vertical and both axes use µm.
    imagesc(ax,x_um,z_um,img.',[0 1]); axis(ax,'xy');
    colormap(ax,flipud(gray(256)));
    xlim(ax,[-5 5]); ylim(ax,[-z_view z_view]);
    daspect(ax,[1 1 1]);
    xticks(ax,[-4 0 4]);
    xlabel(ax,'X (µm)'); ylabel(ax,'Z (µm)');
    title(ax,heading);
end
