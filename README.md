# MATLAB Discrete Trigonometric Transform Library
[![View MATLAB Discrete Trigonometric Transform Library on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/75071-matlab-discrete-trigonometric-transform-library)

## Author

This library is written by Bradley Treeby, University College London. Contact: b.treeby@ucl.ac.uk

## Overview

This repository provides a MATLAB interface to the FFTW implementations of the discrete cosine transforms (DCT) and discrete sine transforms (DST). The functions are implemented in C++ as MATLAB mex functions.

MATLAB natively includes an interface to the complex-to-complex transforms in FFTW via the inbuilt functions `fft`, `fft2`, and `fftn`. However, MATLAB does not include an interface to the real-to-real transforms which correspond to discrete trigonometric transforms (DTTs). This library is intended to fill the gap.

Three functions are included: `dtt1D`, `dtt2D`, and `dtt3D`. These compute DTTs in 1D, 2D, and 3D. The function `dtt1D` can also perform 1D transformations over 2D arrays.

The type of DTT is specified by the input `dtt_type`, where 1 to 4 corresponds to DCTs, and 5 to 8 to DSTs. Currently, only double precisions transforms are supported, thus the input array must be in double precision.

The DTT functions in this library create the FFTW plan using FFTW_ESTIMATE, and then destroy the plan after execution. Consequently, the performance will not match that of the MATLAB inbuilt FFT functions which can re-use wisdom (see `help fftw`).

## Compilation

The mex files can be compiled using `compileDttMex`, which calls `mex`. In many cases, the default paths will need to be changed. Open `compileDttMex` for further details on how to setup and compile. 

Precompiled mex files for Windows and macOS are included in the repository. These have been compiled using MATLAB 2019b. The Windows mex functions were compiled using Windows 10 (1803) and Microsoft Visual C++ 2015. The macOS mex functions were compiled using macOS Catalina (10.15.3) and Xcode Clang++ (11.3.1).

## Examples

An example of using `dtt1D` is included in the function `gradientDTT1D`. This computes a spectral gradient using any of the eight supported DTT symmetries, including an option for grid staggering. Several other example scripts are also included in the examples folder. 

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Change Log

* v1.1 (21 April 2020): 
  * Fixed bug in `gradientDtt1D` (missing `numDim` function)
  * Added `align_output` option to `gradientDtt1D`
  * Updated `example_wave_eq_pstd_1D` to enforce correct boundary position
  * Added `example_wave_eq_pstd_1D_non_reflecting` to simulate general boundary conditions
  * Added `example_wsws_gradient` to illustrate grid shifting with and without `align_output` 
* Added `example_wave_eq_pstd_2D_neumann` and `example_wave_eq_pstd_2D_dirichlet`  
  
* v1.0 (17 April 2020): 
  * Initial release