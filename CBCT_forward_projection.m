%% CBCT forward projection demo (GPU version)
clc; clear; close all;

%% 1) Initialize GPU
g = gpuDevice;
reset(g);

%% 2) Load 3D phantom
filename = '3D_phantom.tif';
Vol = single(imstackread(filename));
Vol = Vol / max(Vol(:));                     % normalize to [0, 1]

[Ny, Nx, Nz] = size(Vol);

%% 3) Define volume coordinates
dx = 1.0; dy = 1.0; dz = 1.0;         % voxel size
step = single(min([dx, dy, dz]));     % ray marching step

x = single(((1:Nx) - (Nx+1)/2) * dx);
y = single(((1:Ny) - (Ny+1)/2) * dy);
z = single(((1:Nz) - (Nz+1)/2) * dz);

%% 4) Define CBCT geometry
SOD = single(200);                    % source-to-object-center distance
SID = single(400);                    % source-to-detector distance

Nu = 256; Nv = 256;                   % detector size
du = single(1.0); dv = single(1.0);   % detector pixel size

u = single(((1:Nu) - (Nu+1)/2) * du);
v = single(((1:Nv) - (Nv+1)/2) * dv);

angles = single(0:0.5:355);           % full scan
Nang = numel(angles);

%% 5) Move fixed data to GPU
Vol_g = gpuArray(Vol); % 3D volume (phantom)
x_g = gpuArray(x); y_g = gpuArray(y); z_g = gpuArray(z); % Coordinates of the volume grid
u_g = gpuArray(u); v_g = gpuArray(v); % Detector coordinates (2D)

proj_gpu = gpuArray.zeros(Nv, Nu, Nang, 'single');

%% 6) Forward projection
for ia = 1:Nang
    theta = single(deg2rad(double(angles(ia))));
    fprintf('Projection %d / %d, angle = %.1f deg\n', ia, Nang, angles(ia));

    proj_gpu(:,:,ia) = cbctForwardOneAngle( ...
        Vol_g, x_g, y_g, z_g, theta, SOD, SID, u_g, v_g, step);
end

wait(g);
proj = gather(proj_gpu);

%% 7) Save results
imstackwrite(uint16(proj * 1000), 'cbct_projections_angle.tif');
save('cbct_projections_angle.mat', 'proj', 'angles');

function P = cbctForwardOneAngle(Vol, x, y, z, theta, SOD, SID, u, v, step)
% Simulate one CBCT projection at angle theta.
Nv = numel(v); Nu = numel(u);

%% 1) Source and detector geometry
src = single([SOD*sin(theta), -SOD*cos(theta), 0]);   % source position
dir0 = -src / norm(src);                              % source -> isocenter
detCenter = src + SID * dir0;                         % detector center

% Detector local axes
ev = single([0, 0, 1]);                               % detector vertical axis
eu = cross(ev, dir0);  eu = eu / norm(eu);            % detector horizontal axis
ev = cross(dir0, eu);  ev = ev / norm(ev);            % re-orthogonalize

%% 2) Detector pixel coordinates
[U, Vv] = meshgrid(u, v);

detX = detCenter(1) + U .* eu(1) + Vv .* ev(1);
detY = detCenter(2) + U .* eu(2) + Vv .* ev(2);
detZ = detCenter(3) + U .* eu(3) + Vv .* ev(3);

%% 3) Ray directions: source -> detector pixels
dirX = detX - src(1);
dirY = detY - src(2);
dirZ = detZ - src(3);

dirNorm = sqrt(dirX.^2 + dirY.^2 + dirZ.^2);
dirX = dirX ./ dirNorm;
dirY = dirY ./ dirNorm;
dirZ = dirZ ./ dirNorm;

%% 4) Ray-box intersection
% Find the common integration range [t_start, t_end]
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
zmin = min(z); zmax = max(z);

epsVal = single(1e-6);
Rx = dirX; Ry = dirY; Rz = dirZ;
Rx(abs(Rx) < epsVal) = epsVal;
Ry(abs(Ry) < epsVal) = epsVal;
Rz(abs(Rz) < epsVal) = epsVal;

tx1 = (xmin - src(1)) ./ Rx;  tx2 = (xmax - src(1)) ./ Rx;
ty1 = (ymin - src(2)) ./ Ry;  ty2 = (ymax - src(2)) ./ Ry;
tz1 = (zmin - src(3)) ./ Rz;  tz2 = (zmax - src(3)) ./ Rz;

tmin = max(max(min(tx1, tx2), min(ty1, ty2)), min(tz1, tz2));
tmax = min(min(max(tx1, tx2), max(ty1, ty2)), max(tz1, tz2));

valid = (tmax > tmin) & (tmax > 0);
P = gpuArray.zeros(Nv, Nu, 'single');

if ~gather(any(valid, 'all'))
    return;
end

t_start = gather(max(single(0), min(tmin(valid), [], 'all')));
t_end   = gather(max(tmax(valid), [], 'all'));

%% 5) Ray integration on GPU
tVec = single(t_start):single(step):single(t_end);

if isempty(tVec)
    return;
end

T = reshape(gpuArray(tVec), 1, 1, []);

Xq = src(1) + T .* dirX;
Yq = src(2) + T .* dirY;
Zq = src(3) + T .* dirZ;

% Trilinear interpolation in the volume
vals = interp3(x, y, z, Vol, Xq, Yq, Zq, 'linear', 0);

% Approximate line integral along each ray
P = sum(vals, 3) * step;
end