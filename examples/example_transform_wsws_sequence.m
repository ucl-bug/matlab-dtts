% DESCRIPTION:
%     This example script defines a symmetric-periodic-sequence with WSWS
%     symmetry, and then uses discrete cosine and sine transforms to
%     compute: 
%
%         1. The gradient
%         2. The gradient interpolated onto a staggered grid
%         3. The function interpolated onto a staggered grid
%
%     This script reproduces the samples used for Figure 4 of [1].
%
%     [1] E. Wise, J. Jaros, B. Cox, and B. Treeby, "Pseudospectral
%     time-domain (PSTD) methods for the wave equation: Realising boundary
%     conditions with discrete sine and cosine transforms", 2020.
%       
% ABOUT:
%     author      - Bradley Treeby
%     date        - 31 March 2020
%     last update - 20 April 2020
%
% See also dtt1D, gradientDtt1D

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

% =========================================================================
% LITERALS
% =========================================================================

% define literals for DTT types used by the function dtt1D 
DCT1                = 1;    % WSWS
DCT2                = 2;    % HSHS
DCT3                = 3;    % WSWA
DCT4                = 4;    % HSHA
DST1                = 5;    % WAWA
DST2                = 6;    % HAHA
DST3                = 7;    % WAWS
DST4                = 8;    % HAHS

% =========================================================================
% INPUT
% =========================================================================

% create sequence with WSWS symmetry
Nx                  = 9;
dx                  = pi/(Nx-1);
x_in                = 0:dx:pi;
f                   = cos(x_in);

% define plot axes
x                   = (0:Nx-1) * dx;
x_sg                = x + dx / 2;

% compute the implied period of the input function
M                   = 2 .* (Nx - 1);

% define DTT wavenumbers for WSWS / DCT-I
n                   = (0:1:M/2);
k                   = 2 .* pi .* n ./ (M .* dx);

% =========================================================================
% GRADIENT
% =========================================================================

% compute forward transform  using DCT-I and multiply by the wavenumbers
f_deriv             = -k .* dtt1D(f, DCT1);
            
% to map from WSWS -> WAWA, remove both endpoints
f_deriv             = f_deriv(2:end - 1);

% compute inverse transform and normalise by the implied period of the
% input function, where S1^-1 = 1/M * S1
f_deriv             = dtt1D(f_deriv, DST1) ./ M;
            
% =========================================================================
% GRADIENT + INTERPOLATION
% =========================================================================
     
% compute forward transform using DCT-I and multiply by the wavenumbers
f_deriv_interp      = -k .* dtt1D(f, DCT1);
            
% to map from WSWS -> HAHA, remove left endpoint
f_deriv_interp      = f_deriv_interp(2:end);

% compute inverse transform and normalise by the implied period of the
% input function, where S2^-1 = 1/M * S3
f_deriv_interp      = dtt1D(f_deriv_interp, DST3) ./ M;

% =========================================================================
% INTERPOLATION
% =========================================================================

% compute forward transform using DCT-I 
f_interp            = dtt1D(f, DCT1);

% to map from WSWS -> HSHS, remove right endpoint
f_interp            = f_interp(1:end - 1);

% compute inverse transform and normalise by the implied period of the
% input function, where C2^-1 = 1/M * C3
f_interp            = dtt1D(f_interp, DCT3) ./ M;

% =========================================================================
% PLOTTING
% =========================================================================

% plot sequence
figure;
subplot(2, 2, 1);
plotData(x, f, 'Sequence (WSWS)');

% plot gradient
subplot(2, 2, 2);
plotData(x(2:end-1), f_deriv, 'Gradient (WAWA)');

% plot gradient + interpolation
subplot(2, 2, 3);
plotData(x_sg(1:end-1), f_deriv_interp, 'Gradient + Interpolation (HAHA)');

% plot interpolation
subplot(2, 2, 4);
plotData(x_sg(1:end-1), f_interp, 'Interpolation (HSHS)');

% nested function to set plot axes, etc, for each subplot
function plotData(x_data, y_data, title_string)
    plot(x_data, y_data, 'k.', 'MarkerSize', 14);
    title(title_string);
    legend(['N = ' num2str(length(y_data))]);
    set(gca, ...
        'XLim', [0, pi], ...
        'YLim', [-1, 1], ...
        'XTick', [0, pi/4, pi/2, 3*pi/4, pi], ...
        'XTickLabel', {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    grid on;
end