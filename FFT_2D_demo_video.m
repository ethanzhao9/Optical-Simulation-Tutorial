%% FFT2 step-by-step demo
clear; clc; close all;

%% ----- Build test image and reference FFT2
N = 10;
x = 0:N-1; y = 0:N-1; [X,Y] = meshgrid(x,y);
fx = 2; fy = 3;
g  = cos(2*pi*fx.*X/N) + 0.6*cos(2*pi*fy.*Y/N);

F2      = fft2(g);
F2s     = fftshift(F2);       % reference 2-D spectrum, centered (both dims)
REFspec = log1p(abs(F2s));    % display-friendly reference

% Centered frequency index for 1-D plots (even N): -N/2 ... N/2-1
k = -N/2 : (N/2-1);

%% ----- Quick context view
figure('Name','Context: image and reference 2-D FFT','Color','w','Position',[80 80 920 410]);
tA = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile(tA,1); imagesc(g); axis image off; colormap gray; colorbar;
title(sprintf('Original 10\\times10  (fx=%d, fy=%d)',fx,fy));
nexttile(tA,2); imagesc(REFspec); axis image off; colormap gray; colorbar;
title('Reference centered 2-D FFT  log(1+|F|)');

%% ----- Video setup
vname = 'fft2_build_from_1d_v4_fixed.mp4';
vw = VideoWriter(vname,'MPEG-4'); vw.FrameRate = 3; vw.Quality = 100; open(vw);

%% ========== ROW PROCESS (build R = FFT_x{g}) ==========
% 8 panels layout (2x4), per-row frame
R = zeros(N,N);          % Row-FFT matrix (unshifted)
F = zeros(N,N);          % Not used for reconstruction in row stage

figR = figure('Name','Row process (8 panels per row)','Color','w','Position',[20 40 1560 820]);
tR = tiledlayout(figR,2,4,'Padding','compact','TileSpacing','compact');

for r = 1:N
    % --- Compute this row's 1-D FFT (along x -> u) and update R ---
    row_signal = g(r,:);               % input (space x) for row r
    row_fft    = fft(row_signal);      % FFT_x{row r}
    row_fft_s  = fftshift(row_fft);    % centered for display
    R(r,:)     = row_fft;              % accumulate into R(y,u)

    % Build a "partial R" that contains only rows 1..r (others zero) for fair reconstruction
    R_partial = zeros(N,N);
    R_partial(1:r, :) = R(1:r, :);

    % --- Current Row-FFT image (center along u ONLY) & its difference to final 2-D FFT ---
    Rvis_centered = fftshift(R, 2);            % center ALONG u ONLY (dim=2)
    Rvis_log      = log1p(abs(Rvis_centered)); % current FFT image (row stage)
    diff_spec_row = log1p(abs(Rvis_centered - F2s));  % diagnostic difference (domains differ by design)

    % --- Row-stage reconstruction: inverse along x only from R_partial ---
    g_rec_row = real(ifft(R_partial, [], 2));  % rows 1..r reconstructed, others ~0
    diff_img  = abs(g_rec_row - g);

    % ----------------- PANELS -----------------

    % [1] Original image — highlight the EXACT row being transformed
    ax1 = nexttile(tR,1); cla(ax1);
    imagesc(ax1, g); axis(ax1,'image'); axis(ax1,'off'); colormap(ax1, gray); colorbar(ax1);
    title(ax1, sprintf('[1] Original image — row %d/%d', r, N));
    hold(ax1,'on'); rectangle(ax1,'Position',[0.5, r-0.5, N, 1], 'EdgeColor','r','LineWidth',1.8); hold(ax1,'off');

    % [2] Input signal (this row, magnitude vs x)
    ax2 = nexttile(tR,2); cla(ax2);
    stem(ax2, 0:N-1, abs(row_signal), 'filled'); grid(ax2,'on');
    xlim(ax2,[-0.5, N-0.5]); xticks(ax2,0:N-1);
    xlabel(ax2,'x (pixel)'); ylabel(ax2,'|row_r(x)|');
    title(ax2, sprintf('[2] Input (row %d) in space', r));

    % [3] FFT of input (centered magnitude vs u)
    ax3 = nexttile(tR,3); cla(ax3);
    stem(ax3, k, abs(row_fft_s), 'filled'); grid(ax3,'on');
    xlim(ax3,[min(k)-0.5, max(k)+0.5]); xticks(ax3,k);
    xlabel(ax3,'u (centered)'); ylabel(ax3,'|FFT_x{row r}|');
    title(ax3, '[3] FFT of input (centered)');

    % [4] Current FFT image (Row-FFT image, centered in u, log)
    ax4 = nexttile(tR,4); cla(ax4);
    imagesc(ax4, Rvis_log); axis(ax4,'image'); axis(ax4,'off'); colormap(ax4, gray); colorbar(ax4);
    title(ax4, '[4] Current FFT image: Row-FFT  log(1+|R|)');

    % [5] Reference 2-D FFT (centered, log)
    ax5 = nexttile(tR,5); cla(ax5);
    imagesc(ax5, REFspec); axis(ax5,'image'); axis(ax5,'off'); colormap(ax5, gray); colorbar(ax5);
    title(ax5, '[5] Reference FFT image  log(1+|F|)');

    % [6] Spectrum difference (Row-FFT vs final 2-D FFT — diagnostic, changes with r)
    ax6 = nexttile(tR,6); cla(ax6);
    imagesc(ax6, diff_spec_row); axis(ax6,'image'); axis(ax6,'off'); colormap(ax6, gray); colorbar(ax6);
    title(ax6, '[6] Spectrum diff  log(1+|Row-FFT - F_ref|)');

    % [7] Reconstruction from current Row info (inverse along x only)
    ax7 = nexttile(tR,7); cla(ax7);
    imagesc(ax7, g_rec_row); axis(ax7,'image'); axis(ax7,'off'); colormap(ax7, gray); colorbar(ax7);
    title(ax7, '[7] Reconstruction from current Row-FFT:  Re{ifft_x(R_partial)}');

    % [8] iFFT image difference |Recon - Original|
    ax8 = nexttile(tR,8); cla(ax8);
    imagesc(ax8, diff_img); axis(ax8,'image'); axis(ax8,'off'); colormap(ax8, gray); colorbar(ax8);
    title(ax8, '[8] |Re{ifft_x(R_partial)} - g|');

    drawnow; frame = getframe(figR); writeVideo(vw, frame);
end


%% ========== COLUMN PROCESS (build F(:,u) = FFT_y{R(:,u}) ==========
% Panel order mirrored with row process (current FFT image in top-right; reference in bottom-left)
figC = figure('Name','Column process (8 panels per column)','Color','w','Position',[40 60 1560 820]);
tC = tiledlayout(figC,2,4,'Padding','compact','TileSpacing','compact');

F(:,:) = 0;                     % start building real 2-D FFT now
R_unshifted = R;
R_shifted   = fftshift(R, 2);   % center ALONG u ONLY (dim=2) for display
Rvis_final  = log1p(abs(R_shifted));

for u_raw = 1:N
    % Map raw column u to shifted display column u_disp (even N):
    u_disp = mod(u_raw-1 + N/2, N) + 1;
    u_centered = (u_raw - 1) - N/2;

    % Column input and its vertical FFT
    col_signal = R_unshifted(:,u_raw);   % input vs y (complex)
    col_fft    = fft(col_signal);        % FFT_y{R(:,u_raw)}
    col_fft_s  = fftshift(col_fft);      % centered
    F(:,u_raw) = col_fft;                % fill current 2-D spectrum (unshifted domain)

    %% Panel [1]: Row-FFT image with active shifted column highlighted
    ax1 = nexttile(tC,1); cla(ax1);
    imagesc(ax1, Rvis_final); axis(ax1,'image'); axis(ax1,'off'); colormap(ax1,gray); colorbar(ax1);
    title(ax1, sprintf('[1] Row-FFT image  (raw u=%d, centered u=%+d)', u_raw, u_centered));
    hold(ax1,'on'); rectangle(ax1,'Position',[u_disp-0.5,0.5,1,N],'EdgeColor','y','LineWidth',2); hold(ax1,'off');

    %% Panel [2]: Input column signal |R(:,u)| vs y
    ax2 = nexttile(tC,2); cla(ax2);
    stem(ax2,0:N-1,abs(col_signal),'filled'); grid(ax2,'on'); xlim(ax2,[-0.5,N-0.5]); xticks(ax2,0:N-1);
    xlabel(ax2,'y (row)'); ylabel(ax2,'|R(y,u)|');
    title(ax2, sprintf('[2] Input (col of R) for u=%d', u_raw));

    %% Panel [3]: FFT of input (centered magnitude vs v)
    ax3 = nexttile(tC,3); cla(ax3);
    stem(ax3,k,abs(col_fft_s),'filled'); grid(ax3,'on'); xlim(ax3,[min(k)-0.5,max(k)+0.5]); xticks(ax3,k);
    xlabel(ax3,'v (centered)'); ylabel(ax3,'|FFT_y{R(:,u)}|');
    title(ax3, '[3] FFT of input (centered)');

    %% Panel [4]: Current FFT image (partial 2-D FFT, centered)
    ax4 = nexttile(tC,4); cla(ax4);
    Fshow = log1p(abs(fftshift(F)));    % center BOTH dims for 2-D spectrum display
    imagesc(ax4,Fshow); axis(ax4,'image'); axis(ax4,'off'); colormap(ax4,gray); colorbar(ax4);
    title(ax4, sprintf('[4] Current FFT image  (cols 1..%d filled)', u_raw));

    %% Panel [5]: Reference FFT image (moved to 2nd row to mirror row layout)
    ax5 = nexttile(tC,5); cla(ax5);
    imagesc(ax5,REFspec); axis(ax5,'image'); axis(ax5,'off'); colormap(ax5,gray); colorbar(ax5);
    title(ax5, '[5] Reference FFT image  log(1+|F|)');

    %% Panel [6]: Spectrum difference (partial vs reference)
    ax6 = nexttile(tC,6); cla(ax6);
    diff_spec = log1p(abs(fftshift(F) - F2s));  % both centered for fair visual comparison
    imagesc(ax6,diff_spec); axis(ax6,'image'); axis(ax6,'off'); colormap(ax6,gray); colorbar(ax6);
    title(ax6, '[6] Spectrum difference  log(1+|F_current - F_ref|)');

    %% Panel [7]: Reconstruction from current F
    ax7 = nexttile(tC,7); cla(ax7);
    g_rec = real(ifft2(F));
    imagesc(ax7,g_rec); axis(ax7,'image'); axis(ax7,'off'); colormap(ax7,gray); colorbar(ax7);
    title(ax7, '[7] Reconstruction from current F:  Re{ifft2(F)}');

    %% Panel [8]: iFFT image difference |Recon - Original|
    ax8 = nexttile(tC,8); cla(ax8);
    imagesc(ax8,abs(g_rec - g)); axis(ax8,'image'); axis(ax8,'off'); colormap(ax8,gray); colorbar(ax8);
    title(ax8, '[8] |Re{ifft2(F)} - g|');

    drawnow; frame = getframe(figC); writeVideo(vw,frame);
end

%% ----- Wrap up
close(vw);
fprintf('Video written to: %s\n', vname);

% Sanity checks
fprintf('||F (built) - fft2(g)||_2 = %.3e\n', norm(F(:)-F2(:)));
