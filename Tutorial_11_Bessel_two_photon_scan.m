% Simulate one XY scan with a conventional or an extended Bessel two-photon focus.
clear; close all; clc;

%% ===================== Parameters =====================
lambda = 920e-9;                   % Vacuum excitation wavelength (m)
n = 1.33;                          % Sample refractive index
NA_obj = 0.65;                     % Shared objective cutoff
NA_gauss = 0.34;                   % 1/e amplitude radius of the Gaussian pupil
NA_ring = 0.32;                    % Annulus mean radius, in NA units
dNA_ring = 0.045;                   % Full radial width of the hard annulus
dxy = 0.25e-6; dz = 0.50e-6;        % Lateral and axial voxel sizes (m)
x = (-216:216)*dxy; y = x;
z = (-100:100)*dz;
n_random = 1000;
r_kernel = 30e-6;                   % Half-width of the PSF used in the scan (m)
output_dir = fullfile(fileparts(mfilename('fullpath')),'two_photon_figures');
if ~exist(output_dir,'dir'), mkdir(output_dir); end

%% ===================== 3-D Sample =====================
% Array order is [Y,X,Z]. Each bead is a unit-brightness point emitter.
rng(0);
ix_pool = find(abs(x*1e6)<=47);
iy_pool = find(y*1e6>=-33 & y*1e6<=47);
iz_pool = find(abs(z*1e6)<=42);
sample_size = [numel(y),numel(x),numel(z)];
random_ids = randperm(numel(iy_pool)*numel(ix_pool)*numel(iz_pool),n_random);
[iy,ix,iz_random] = ind2sub([numel(iy_pool),numel(ix_pool),numel(iz_pool)],random_ids);
bead_ids = sub2ind(sample_size,iy_pool(iy),ix_pool(ix),iz_pool(iz_random));

% The last two beads share XY but lie on opposite sides of focus.
control_xyz = [-32 -43 0; -16 -43 14; 0 -43 27; 16 -43 40; 32 -43 -14; 32 -43 14]; % um
control_idx = round(control_xyz*1e-6./[dxy dxy dz])+([numel(x) numel(y) numel(z)]+1)/2;
control_ids = sub2ind(sample_size,control_idx(:,2),control_idx(:,1),control_idx(:,3));
object = false(sample_size);
object([bead_ids(:);control_ids]) = true;
[by,bx,bz] = ind2sub(sample_size,find(object));
iz0 = (numel(z)+1)/2;
assert(nnz(object)==n_random+size(control_xyz,1),'Bead count mismatch.');

%% ===================== Two-Photon Excitation PSFs =====================
% Circularly symmetric pupils: scalar angular spectrum integrated over NA.
% For two-photon excitation, H = intensity^2 = abs(field)^4.
half_kernel = round(r_kernel/dxy);
q = (-half_kernel:half_kernel)*dxy;
[Xq,Yq] = meshgrid(q,q);
rq = hypot(Xq,Yq);
r = (0:ceil(max(rq(:))/(dxy/8)))*(dxy/8);
a_gauss = linspace(0,NA_obj,601);
p_gauss = exp(-(a_gauss/NA_gauss).^2);
a_ring = linspace(NA_ring-dNA_ring/2,NA_ring+dNA_ring/2,301);
assert(a_ring(1)>0 && a_ring(end)<=NA_obj,'Annulus must fit the pupil.');
Hg_r = two_photon_psf(r,z,lambda,n,a_gauss,p_gauss);
Hb_r = two_photon_psf(r,z,lambda,n,a_ring,ones(size(a_ring)));

% Build [Y,X,Z] kernels on the same grid as the sample.
Hg = reshape(single(interp1(r(:),Hg_r,rq(:),'linear')),[size(rq),numel(z)]);
Hb = reshape(single(interp1(r(:),Hb_r,rq(:),'linear')),[size(rq),numel(z)]);
center = half_kernel+1;
g_xy = width_at_half(r,Hg_r(:,iz0));
b_xy = width_at_half(r,Hb_r(:,iz0));
g_z = width_at_half(z(iz0:end),Hg_r(1,iz0:end));
b_z = width_at_half(z(iz0:end),Hb_r(1,iz0:end));


%% ===================== Image Formation at Z = 0 =====================
G = zeros(numel(y),numel(x),'single');     % Conventional scan
B = zeros(size(G),'single');              % Bessel scan
for iz = 1:numel(z)
    sample = single(object(:,:,iz));    % Beads at z(iz)
    j = numel(z)-iz+1;                   % PSF at 0 - z(iz)
    G = G + conv2(sample,Hg(:,:,j),'same');
    B = B + conv2(sample,Hb(:,:,j),'same');
end

%% ===================== Read the Result =====================
control_depths = control_idx(:,3);
control_g = Hg_r(1,control_depths);
control_b = Hb_r(1,control_depths);
fprintf('Beads: %d random + %d controls. Raster: %d x %d pixels, one scan each.\n',n_random,size(control_xyz,1),numel(x),numel(y));
fprintf('Gaussian/Bessel FWHM XY: %.3f / %.3f um; Z: %.2f / %.2f um.\n',g_xy*1e6,b_xy*1e6,g_z*1e6,b_z*1e6);
random_depth = z(iz_pool(iz_random));
fprintf('Random bead centers in the on-axis half-maximum interval: %d Gaussian, %d Bessel.\n', ...
    nnz(abs(random_depth)<=g_z/2),nnz(abs(random_depth)<=b_z/2));
fprintf('Control response (depth um, Gaussian, Bessel):\n');
disp([control_xyz(:,3),control_g(:),control_b(:)]);
fprintf('Isolated pair, predicted on-axis sum: Gaussian %.5f; Bessel %.5f.\n', ...
    sum(control_g(5:6)),sum(control_b(5:6)));
image_max = double(max([G(:);B(:)]));
fprintf('Shared image peak: %.3f (single in-focus bead peak = 1).\n',image_max);
um = 1e6;

% Figure 1: the sample truth and the two measured rasters.
f = figure('Color','w','Position',[50 50 1540 560]);
t = tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
a = nexttile(t);
scatter(a,x(bx)*um,y(by)*um,11,z(bz)*um,'filled','MarkerFaceAlpha',0.84);
axis(a,'image'); xlim(a,[-54 54]); ylim(a,[-54 54]); box(a,'on');
xlabel(a,'X (µm)'); ylabel(a,'Y (µm)'); title(a,'A   Bead locations (XY)');
colormap(a,turbo(256)); caxis(a,[-42 42]);
cb = colorbar(a,'southoutside'); cb.Label.String = 'Bead depth Z (µm)';
for k = 1:2
    a = nexttile(t);
    if k==1
        xy_image(a,x*um,y*um,G,image_max,'B   Conventional two-photon');
    else
        xy_image(a,x*um,y*um,B,image_max,'C   Bessel two-photon');
    end
    cb = colorbar(a,'southoutside'); cb.Label.String = 'Relative fluorescence';
end
title(t,'Same 3D sample, same XY raster at Z = 0','FontWeight','bold');
exportgraphics(f,fullfile(output_dir,'01_same_sample_two_scans.png'),'Resolution',180);

% Figure 2: the axial excitation extent (same physical aspect ratio).
f = figure('Color','w','Position',[50 50 1080 555]);
t = tiledlayout(f,1,4,'Padding','compact','TileSpacing','compact');
sel = abs(q*um)<=13;
for k = 1:2
    a = nexttile(t);
    if k==1
        map = squeeze(Hg(center,sel,:)).'; heading = 'A   Conventional focus';
    else
        map = squeeze(Hb(center,sel,:)).'; heading = 'B   Bessel focus';
    end
    imagesc(a,q(sel)*um,z*um,map,[0 1]);
    axis(a,'image'); axis(a,'xy'); colormap(a,hot(256));
    xlabel(a,'X (µm)'); ylabel(a,'Z (µm)'); title(a,heading);
    cb = colorbar(a); cb.Label.String = 'Two-photon excitation';
end
a = nexttile(t,[1 2]); hold(a,'on');
plot(a,z*um,Hg_r(1,:),'LineWidth',2,'Color',[0.11 0.35 0.72]);
plot(a,z*um,Hb_r(1,:),'LineWidth',2,'Color',[0.8 0.27 0.06]);
yline(a,0.5,':','HandleVisibility','off'); box(a,'on');
xlabel(a,'Z (µm)'); ylabel(a,'On-axis response'); ylim(a,[0 1.05]);
xlim(a,[z(1) z(end)]*um); title(a,'C   Along the optical axis');
legend(a,{sprintf('Conventional: %.1f µm',g_z*um),sprintf('Bessel: %.1f µm',b_z*um)}, ...
    'Location','southoutside');
title(t,'Two-photon excitation PSF | Y = 0','FontWeight','bold');
exportgraphics(f,fullfile(output_dir,'02_excitation_psfs.png'),'Resolution',180);

% Figure 3: enlarge the same control strip from both simulated scans.
f = figure('Color','w','Position',[50 50 1120 355]);
t = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
for k = 1:2
    a = nexttile(t);
    if k==1
        xy_image(a,x*um,y*um,G,image_max,'Conventional two-photon');
    else
        xy_image(a,x*um,y*um,B,image_max,'Bessel two-photon');
    end
    xlim(a,[-40 40]); ylim(a,[-48 -36]); xticks(a,control_xyz(1:5,1));
    if k==1
        a.XTickLabel = []; xlabel(a,'');
        for j = 1:5
            if j==1
                label = 'Z = 0 µm';
            elseif j<5
                label = sprintf('Z = +%d µm',control_xyz(j,3));
            else
                label = 'Z = ±14 µm';
            end
            text(a,control_xyz(j,1),-38,label,'Color','w', ...
                'HorizontalAlignment','center','FontSize',11);
        end
    end
end
title(t,'Control beads | the rightmost two share the same XY position','FontWeight','bold');
exportgraphics(f,fullfile(output_dir,'03_control_beads.png'),'Resolution',180);

%% ===================== 3-D Scan Coverage =====================
% Sweep the on-axis half-maximum segment across the XY raster. The resulting
% slab marks a response threshold; it is not a mask used in the simulation.
half_bessel_um = b_z*um/2;
fprintf('Bessel scan coverage: X [%.2f, %.2f], Y [%.2f, %.2f], Z [%.2f, %.2f] um.\n', ...
    x([1 end])*um,y([1 end])*um,[-half_bessel_um half_bessel_um]);
f = figure('Color','w','Position',[50 50 1170 785]);
a = axes(f,'Position',[0.08 0.15 0.84 0.74]); hold(a,'on');
[yr,xr,zr] = ind2sub(sample_size,bead_ids);
x_random = x(xr)*um; y_random = y(yr)*um; z_random = z(zr)*um;
in_band = abs(z_random)<=half_bessel_um;
slab_color = [0.08 0.45 0.67];
pair_color = [0.58 0.12 0.20];
beam_color = [0.85 0.36 0.04];
x_lim = x([1 end])*um; y_lim = y([1 end])*um;
rim_x = x_lim([1 2 2 1 1]); rim_y = y_lim([1 1 2 2 1]);
for z_face = [-half_bessel_um half_bessel_um]
    patch(a,rim_x,rim_y,z_face*ones(1,5),slab_color,'FaceAlpha',0.08, ...
        'EdgeColor',slab_color,'LineStyle','--','LineWidth',1.3,'HandleVisibility','off');
end
for k = 1:4
    plot3(a,[rim_x(k) rim_x(k)],[rim_y(k) rim_y(k)],[-half_bessel_um half_bessel_um], ...
        '--','Color',slab_color,'LineWidth',1.3,'HandleVisibility','off');
end
focus = plot3(a,rim_x,rim_y,zeros(1,5),'-','Color',[0.24 0.27 0.31],'LineWidth',1.2);
outside = scatter3(a,x_random(~in_band),y_random(~in_band),z_random(~in_band),15, ...
    [0.48 0.52 0.58],'filled','MarkerFaceAlpha',0.58);
inside = scatter3(a,x_random(in_band),y_random(in_band),z_random(in_band),19, ...
    slab_color,'filled','MarkerFaceAlpha',0.83);
% Use the exact scan kernel at the XY address shared by the control pair.
% The 50% isosurface shows the central lobe on the same physical axes.
beam_xy_um = control_xyz(5,1:2);
beam_roi = abs(q*um)<=2;
beam_mesh = isosurface(q(beam_roi)*um,q(beam_roi)*um,z*um,Hb(beam_roi,beam_roi,:),0.5);
beam_mesh_world = beam_mesh;
beam_mesh_world.vertices = beam_mesh.vertices + [beam_xy_um 0];
beam = patch(a,beam_mesh_world,'FaceColor',beam_color,'EdgeColor','none');
control = scatter3(a,control_xyz(:,1),control_xyz(:,2),control_xyz(:,3),57,'d', ...
    'MarkerEdgeColor',pair_color,'LineWidth',1.65);
for z_face = [-half_bessel_um half_bessel_um]
    text(a,x_lim(2)+2,y_lim(1),z_face,sprintf('%+.1f µm',z_face), ...
        'Color',slab_color,'FontSize',13,'VerticalAlignment','middle','Clipping','off', ...
        'BackgroundColor','w','Margin',2);
end
xlim(a,x_lim+[-4 12]); ylim(a,y_lim+[-4 4]); zlim(a,z([1 end])*um);
daspect(a,[1 1 1]); view(a,[-34 23]); grid(a,'on'); box(a,'on');
set(a,'FontSize',13,'TickDir','out','GridAlpha',0.13,'Projection','orthographic');
xlabel(a,'X (µm)'); ylabel(a,'Y (µm)'); zlabel(a,'Z (µm)');
title(a,'Bessel focus and scan coverage','FontSize',18);
legend(a,[inside outside focus control beam], ...
    {sprintf('Within half-maximum range (%d)',nnz(in_band)), ...
    sprintf('Outside (%d)',nnz(~in_band)),'XY raster at Z = 0','Control beads', ...
    'Single-focus PSF (50%)'}, ...
    'Location','southoutside','Orientation','horizontal','FontSize',10,'Box','off');

exportgraphics(f,fullfile(output_dir,'04_bessel_scan_volume.png'),'Resolution',180);

% Same 16-bit conversion for both images, with no per-image normalization.
imwrite(uint16(double(G)*(65535/image_max)),fullfile(output_dir,'conventional_xy.tif'));
imwrite(uint16(double(B)*(65535/image_max)),fullfile(output_dir,'bessel_xy.tif'));

%% ===================== Local Functions =====================
function H = two_photon_psf(r,z,lambda,n,a,pupil)
    % Radial area element a*da, with a = n*sin(theta) = lambda*f_r.
    weight = pupil.*a;
    weight([1 end]) = weight([1 end])/2;
    kz = 2*pi/lambda*sqrt(n^2-a.^2);
    J = besselj(0,(2*pi/lambda)*r(:)*a);
    phase = exp(1i*(kz(:)-kz(1)).*z);
    E = (J.*weight)*phase/sum(weight);
    I = abs(E).^2; H = I.^2;
end

function w = width_at_half(q,p)
    i = find(p<=p(1)/2,1);
    assert(~isempty(i) && i>1,'PSF never reached half maximum on the grid.');
    w = 2*interp1(p(i-1:i),q(i-1:i),p(1)/2);
end

function xy_image(a,x,y,img,peak,heading)
    imagesc(a,x,y,img,[0 peak]);
    % A shared brightened display lookup table; TIFF data remain linear.
    axis(a,'image'); axis(a,'xy'); colormap(a,gray(256).^0.65);
    xlim(a,[-54 54]); ylim(a,[-54 54]);
    xlabel(a,'X (µm)'); ylabel(a,'Y (µm)'); title(a,heading);
end
