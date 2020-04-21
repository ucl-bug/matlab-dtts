% DESCRIPTION:
%     This example script shows a WSWS gradient calculated on regular and
%     staggered grids, both with and without alignment.
%
%     Further details are given in [1].
%
%     [1] E. Wise, J. Jaros, B. Cox, and B. Treeby, "Pseudospectral
%     time-domain (PSTD) methods for the wave equation: Realising boundary
%     conditions with discrete sine and cosine transforms", 2020.
% 
% ABOUT:
%     author      - Bradley Treeby
%     date        - 20 April 2020
%     last update - 20 April 2020
%
% Copyright (C) 2020-2020 Bradley Treeby
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

clearvars;

% create function with WSWS symmetry
Nx        = 10;
dx        = pi/(Nx-1);
x_in      = 0:dx:pi;
f         = cos(x_in);

% =========================================================================
% GRADIENT WITH ALIGNMENT
% =========================================================================

% define x vectors
x_out_def = x_in;
x_out_pos = x_in + dx/2;
x_out_neg = x_in - dx/2;

% compute gradient using DCT-I (WSWS symmetry)
dfdx_def  = gradientDtt1D(f, dx, 1, 0);
dfdx_pos  = gradientDtt1D(f, dx, 1, 1);
dfdx_neg  = gradientDtt1D(f, dx, 1, 2);

% analytical derivative
x_an      = -dx:dx/5:2*pi - dx;
dfdx_an   = -sin(x_an);

% plot and compare
figure;
plot(...
    x_an, dfdx_an, 'k-', ...
    x_out_def, dfdx_def, 'r.', ...
    x_out_pos, dfdx_pos, 'bo', ...
    x_out_neg, dfdx_neg, 'g.');
legend(...
    'Analytical Derivative', ...
    'DTT shift = 0', ...
    'DTT shift = 1', ...
    'DTT shift = 2', ...
    'Location', 'Best');
xlabel('x');
ylabel('f''(x)');
title('DCT-I WSWS (align output = true)');
set(gca, 'XLim', [-dx, 2*pi - dx]);
grid on;

% =========================================================================
% GRADIENT WITHOUT ALIGNMENT
% =========================================================================

% define x vectors
x_out_def = x_in(2:end-1);
x_out_pos = x_in(1:end-1) + dx/2;
x_out_neg = x_in(1:end-1) + dx/2;

% compute gradient using DCT-I (WSWS symmetry)
dfdx_def  = gradientDtt1D(f, dx, 1, 0, false);
dfdx_pos  = gradientDtt1D(f, dx, 1, 1, false);
dfdx_neg  = gradientDtt1D(f, dx, 1, 2, false);

% analytical derivative
x_an      = -dx:dx/5:2*pi - dx;
dfdx_an   = -sin(x_an);

% plot and compare
figure;
plot(...
    x_an, dfdx_an, 'k-', ...
    x_out_def, dfdx_def, 'r.', ...
    x_out_pos, dfdx_pos, 'bo', ...
    x_out_neg, dfdx_neg, 'g.');
legend(...
    'Analytical Derivative', ...
    'DTT shift = 0', ...
    'DTT shift = 1', ...
    'DTT shift = 2', ...
    'Location', 'Best');
xlabel('x');
ylabel('f''(x)');
title('DCT-I WSWS (align output = false)');
set(gca, 'XLim', [-dx, 2*pi - dx]);
grid on;