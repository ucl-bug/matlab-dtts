%DTT1D Discrete trigonometric transform.
%
% DESCRIPTION:
%     dtt1D computes the one-dimensional discrete trigonometric transform
%     (DTT) of the input array x using FFTW (http://www.fftw.org). 
%
%     MATLAB natively includes an interface to the complex-to-complex
%     transforms in FFTW via the inbuilt functions fft, fft2, and fftn.
%     However, MATLAB does not include an interface to the real-to-real
%     transforms which correspond to  discrete trigonometric transforms,
%     i.e., discrete cosine transforms (DCT) and discrete sine transforms
%     (DST). This library is intended to fill the gap.
%
%     The type of DTT is specified by the input dtt_type, where 1 to 4
%     corresponds to DCTs, and 5 to 8 to DSTs. Currently, only double
%     precisions transforms are supported, thus the input array must be in
%     double precision.
%
%     The DTT functions in this library create the FFTW plan using
%     FFTW_ESTIMATE, and then destroy the plan after execution.
%     Consequently, the performance will not match that of the inbuilt fft
%     functions which can re-use wisdom (see help fftw).
%
%     For compilation instructions, see compileDttMex.
%
% USAGE:
%     X = dtt1D(x, dtt_type)
%     X = dtt1D(x, dtt_type, dim)
%
% INPUTS:
%     x             - 1D or 2D array to transform in double precision.
%     dtt_type      - Type of discrete trigonometric transform. This
%                     relates to the assumed symmetry of the input x, where
%                     the symmetry can be either whole (W) or half (H)
%                     sample symmetric (S) or antisymmetric (A) at each end
%                     of sample. The dtt_type is specified as an integer
%                     between 1 and 8, corresponding to the 8 DTTs
%                     implemented in FFTW.
%
%                         1: DCT-I    WSWS
%                         2: DCT-II   HSHS
%                         3: DCT-III  WSWA
%                         4: DCT-IV   HSHA
%                         5: DST-I    WAWA
%                         6: DST-II   HAHA
%                         7: DST-III  WAWS
%                         8: DST-IV   HAHS
%
%                     The first four transforms correspond to discrete
%                     cosine transforms, and the second four transforms to
%                     discrete sine transforms.
%     dim           - For 2D arrays, dim specifies the dimension over which
%                     the transform is taken (equivalent to the dim in put
%                     in MATLAB's fft).
%
% OUTPUTS:
%     X             - Discrete trigonometric transform of the input array
%                     x.
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 31 May 2012
%     last update   - 19 March 2020
%
% Copyright (C) 2012-2020 Bradley Treeby
%
% See also dtt2D, dtt3D

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program.  If not, see <https://www.gnu.org/licenses/>.