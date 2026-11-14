# MATLAB Optical Simulation Tutorials

MATLAB examples for Fourier optics, microscopy, and computed tomography. Follow the numbered tutorials to go from a pupil function to simulated microscope images, or use the index to find a specific technique.

## Quick start

1. Download the repository from **Code → Download ZIP**, or clone it with Git.
2. Open the repository folder in MATLAB. Start with [Tutorial 1](Tutorial_1_Simulate_2D_PSF.m) and click **Run** to calculate a focal-plane point spread function (PSF) from a circular pupil.
3. Change the parameters near the top of the script and rerun it. Tutorials 2–4 extend the PSF model to axial propagation, confocal detection, and image formation.

Use MATLAB R2021b or later as a starting point. Extra toolbox, data, and hardware requirements are listed below.

## Tutorial index

| No. | MATLAB script | What it shows | Zhihu article (中文) |
| --- | --- | --- | --- |
| **1–4** | **PSF and confocal imaging** | | |
| 1 | [2D PSF](Tutorial_1_Simulate_2D_PSF.m) | Compute the in-focus intensity from a circular pupil. | [理解并模拟 2D PSF](https://zhuanlan.zhihu.com/p/1902090402686825172) |
| 2 | [3D PSF](Tutorial_2_Simulate_3D_PSF.m) | Propagate the field through focus and inspect lateral and axial sections. | [模拟 3D PSF](https://zhuanlan.zhihu.com/p/1902293550323957913) |
| 3 | [Confocal PSF](Tutorial_3_Simulate_confocal_PSF.m) | Combine excitation and detection PSFs and compare the effect of a pinhole. | [宽场和共聚焦显微镜的 PSF](https://zhuanlan.zhihu.com/p/1902427509158581857) |
| 4 | [Image formation with a pinhole](Tutorial_4_Simulate_PSF_with_pinhole.m) | Image a simulated 3D object with and without a pinhole, including Poisson noise. | [利用 PSF 模拟显微镜成像](https://zhuanlan.zhihu.com/p/1903927228938560123) |
| **5, 8** | **Structured illumination (SIM)** | | |
| 5 | [Moiré fringes](Tutorial_5_Simulate_Moire_Fringes.m) | Visualize how patterned illumination shifts sample frequencies. | [SIM 成像原理](https://zhuanlan.zhihu.com/p/1911024204633071916) |
| 8 | [Raw SIM images](Tutorial_8_SIM_raw_images_simulation.m) | Generate raw frames at multiple illumination phases and orientations and export TIFF stacks. | [SIM 超分辨成像模拟](https://zhuanlan.zhihu.com/p/1932564385773121657) |
| **6–7** | **Phase retrieval and wavefront aberrations** | | |
| 6 | [Phase retrieval](Tutorial_6_Phase_retrieval_demo.m) | Design a phase-only hologram for a 3 × 3 spot array using Gerchberg–Saxton iterations. | [Phase Retrieval 相位恢复](https://zhuanlan.zhihu.com/p/1925315376465543692) |
| 7 | [Zernike aberrations](Tutorial_7_Wavefront_Zernike_function.m) | Apply Zernike phase terms in the pupil and compare the resulting focal PSFs. | [Zernike 多项式模拟波前畸变](https://zhuanlan.zhihu.com/p/1927307947966436540) |
| **9–10** | **STED microscopy** | | |
| 9 | [STED doughnut beam](Tutorial_9_STED_beam_demo.m) | Compare analytical doughnut intensity profiles for different beam waists. | [STED 超分辨成像原理](https://zhuanlan.zhihu.com/p/1942540299667740543) |
| 10 | [2D STED PSF](Tutorial_10_STED_PSF_2D.m) | Simulate a vortex-generated depletion beam and compare effective PSFs at several depletion strengths. | [2D STED 超分辨成像模拟](https://zhuanlan.zhihu.com/p/1981838353243596699) |
| **11** | **Bessel beams and two-photon imaging** | | |
| 11a | [3D Bessel beam](Tutorial_11_Bessel_beam_3D.m) | Propagate a hard annular pupil at 920 nm and compare ring radius and width in XY and XZ. | [贝塞尔光模拟](https://zhuanlan.zhihu.com/p/2104829994572816626) |
| 11b | [Bessel two-photon scan](Tutorial_11_Bessel_two_photon_scan.m) | Scan the same 3D bead sample once in XY with a conventional focus and once with a Bessel focus. | [贝塞尔光模拟（双光子扫描）](https://zhuanlan.zhihu.com/p/2104829994572816626) |

## Additional demos

| MATLAB script | What it shows |
| --- | --- |
| [2D FFT animation](FFT_2D_demo_video.m) | Builds a 2D Fourier transform from 1D row and column transforms and exports an MP4. |
| [Unfiltered backprojection](CT_Unfiltered_Backprojection.m) | Reconstructs a 2D phantom from a sinogram and compares the manual implementation with MATLAB's built-in result. |
| [Filtered backprojection](CT_Filtered_Backprojection.m) | Compares unfiltered backprojection with Ram–Lak and Shepp–Logan filtering. |
| [Cone-beam CT forward projection](CBCT_forward_projection.m) | Simulates cone-beam projections of a 3D TIFF phantom on a GPU. See the requirements below before running. |

## Requirements

- **Tutorial 4** uses Statistics and Machine Learning Toolbox for Poisson noise.
- **Tutorial 6 and the 2D CT demos** use Image Processing Toolbox.
- **Cone-beam CT** requires Parallel Computing Toolbox, a supported GPU, a multipage TIFF named `3D_phantom.tif`, and `imstackread` and `imstackwrite` on the MATLAB path. The phantom and those helper functions need to be supplied separately.

## Contact

For questions about a script, [open an issue](https://github.com/ethanzhao9/Optical-Simulation-Tutorial/issues) or email [Ethan Zhao](mailto:zhaoyc9@163.com).
