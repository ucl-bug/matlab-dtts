% DESCRIPTION:
%     This example script solves the 2D wave equation (written as two
%     coupled first-order equations) using a DTT-based PSTD method subject
%     to Neumann boundary conditions on each side of the domain. 
%
%     The spatial gradients are calculated using discrete cosine transforms
%     (DCTs) and discrete sine transforms (DST) chosen to enfore the
%     required boundary condition. Time integration is performed using a
%     first-order accurate foward difference scheme.
%
%     Due to the symmetry of the DTTs used (WSWS), the position of the
%     boundary is at the first and last grid points in the domain.
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
%     last update - 21 April 2020
%
% Copyright (C) 2020 Bradley Treeby
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

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% transform types used by dtt1D
DCT1 = 1;    % WSWS
DCT2 = 2;    % HSHS
DCT3 = 3;    % WSWA
DCT4 = 4;    % HSHA
DST1 = 5;    % WAWA
DST2 = 6;    % HAHA
DST3 = 7;    % WAWS
DST4 = 8;    % HAHS

% plot frequency for snapshots of the field (time steps)
plot_freq = 250;

% =========================================================================
% DEFINE SIMULATION SETTINGS
% =========================================================================

% set the grid size (assuming grid is square)
Nx   = 256;     % grid size [m]
dx   = 1/Nx;    % grid spacing [m]

% set the medium properties
c0   = 1500;    % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% set CFL and number of time steps
CFL  = 0.2;
Nt   = 1000;

% =========================================================================
% DEFINE INITIAL CONDITIONS
% =========================================================================

% spatial grid
x = (0:Nx - 1) * dx;

% properties of Gaussian initial condition
width = Nx * dx / 28;
offset = Nx * dx / 3;

% define initial pressure distribution on regular grid as a Gaussian
p = 0.5 * exp(-((x - offset) / width).^2);
p = p.' * p;

% define initial velocity to be zero
ux = 0;
uy = 0;

% =========================================================================
% DEFINE DTT VARIABLES
% =========================================================================
        
% compute the implied period of the input function
M = 2 .* (Nx - 1);

% define wavenumber vectors for WSWS / DCT-I
k_wsws = 2 .* pi .* (0:1:M/2).' ./ (M .* dx);

% define wavenumber vectors for HAHA / DST-II
k_haha = 2 .* pi .* (1:1:M/2).' ./ (M .* dx);
      
% =========================================================================
% RUN SIMULATIONS USING DTT-BASED PSTD METHOD
% =========================================================================

% calculate the time step
dt = CFL * dx / c0;

% create custom colour map
cm = [bone(256); flipud(hot(256))];

% calculate pressure in a loop
for time_ind = 1:Nt

    % ---------------
    % calculate dpdx
    % ---------------
    
    % forward transform over x
    dpdx = dtt1D(p, DCT1, 1);
    
    % spectral derivative
    dpdx = bsxfun(@times, -k_wsws, dpdx);
    
    % WSWS -> HAHA, remove left endpoint
    dpdx = dpdx(2:end, :);
    
    % inverse transform over x, S2^-1 = S3
    dpdx = dtt1D(dpdx, DST3, 1) ./ M;
    
    % ---------------
    % calculate dpdy
    % ---------------
    
    % forward transform over y
    dpdy = dtt1D(p, DCT1, 2);
    
    % spectral derivative
    dpdy = bsxfun(@times, -k_wsws.', dpdy);
    
    % WSWS -> HAHA, remove left endpoint
    dpdy = dpdy(:, 2:end);
    
    % inverse transform over y, S2^-1 = S3
    dpdy = dtt1D(dpdy, DST3, 2) ./ M;
    
    % ---------------
    % update u
    % ---------------
    
    % update components of particle velocity
    ux = ux - dt / rho0 * dpdx;
    uy = uy - dt / rho0 * dpdy;

    % ---------------
    % calculate duxdx
    % ---------------
    
    % forward transform over x
    duxdx = dtt1D(ux, DST2, 1);
    
    % spectral derivative
    duxdx = bsxfun(@times, k_haha, duxdx);
    
    % HAHA -> WSWS, prepend zero
    duxdx = [zeros(1, Nx); duxdx]; %#ok<AGROW>
    
    % inverse transform over x, C1^-1 = C1
    duxdx = dtt1D(duxdx, DCT1, 1) ./ M;
    
    % ---------------
    % calculate duydy
    % ---------------
    
    % forward transform over y
    duydy = dtt1D(uy, DST2, 2);
    
    % spectral derivative
    duydy = bsxfun(@times, k_haha.', duydy);
    
    % HAHA -> WSWS, prepend zero
    duydy = [zeros(Nx, 1), duydy]; %#ok<AGROW>
    
    % inverse transform over y, C1^-1 = C1
    duydy = dtt1D(duydy, DCT1, 2) ./ M;      
    
    % ---------------
    % update p
    % ---------------
    
    % update p
    p = p - dt * rho0 * c0^2 * (duxdx + duydy);

    % ---------------
    % plot
    % ---------------    
    
    % plot snapshots of pressure field
    if (time_ind == 1) || ~rem(time_ind, plot_freq)
        figure;
        imagesc(x, x, p, 0.075 * [-1, 1]);
        axis image;
        colormap(cm);
        box on;
        set(gca, 'XTick', 0:0.25:1, 'YTick', 0:0.25:1, 'XTickLabel', {}, 'YTickLabel', {});
        drawnow;
    end

end