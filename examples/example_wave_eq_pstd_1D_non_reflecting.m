% DESCRIPTION:
%     This example script solves the 1D wave equation (written as two
%     coupled first-order equations) using a DTT-based PSTD method subject
%     to four different combinations of Dirichlet and Neumann boundary
%     conditions at each end of the domain. 
%
%     The solutions for different boundary conditions are then summed to
%     give different reflection coefficients at each end of the domain,
%     including non-reflecting boundaries. This builds on
%     example_wave_eq_pstd.m. 
%
%     Further details are given in [1].
%
%     [1] E. Wise, J. Jaros, B. Cox, and B. Treeby, "Pseudospectral
%     time-domain (PSTD) methods for the wave equation: Realising boundary
%     conditions with discrete sine and cosine transforms", 2020.
%       
% ABOUT:
%     author      - Bradley Treeby
%     date        - 19 April 2020
%     last update - 20 April 2020
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

%% ========================================================================
% DEFINE LITERALS
% =========================================================================

% shift and align variables used by gradientDtt1D (shift = 1 and shift = 2
% return the same output if align_output = false)
shift_output = 1;
align_output = false;

% transform types used by gradientDtt1D
DCT1 = 1;    % WSWS
DCT2 = 2;    % HSHS
DCT3 = 3;    % WSWA
DCT4 = 4;    % HSHA
DST1 = 5;    % WAWA
DST2 = 6;    % HAHA
DST3 = 7;    % WAWS
DST4 = 8;    % HAHS

% number of time snapshots to include in waterfall plot
waterfall_snapshots = 80;

% number of boundary conditions
bc_num = 4;

%% ========================================================================
% DEFINE SIMULATION SETTINGS
% =========================================================================

% set the grid size
Nx   = 256;     % grid size [m]
dx   = 1/Nx;    % grid spacing [m]

% set the medium properties
c0   = 1500;    % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% set CFL and number of time steps
CFL  = 0.2;
Nt   = 2720;

%% ========================================================================
% DEFINE INITIAL CONDITIONS
% =========================================================================

% spatial grid
x = (0:Nx - 1) * dx;

% properties of Gaussian initial condition
width = Nx * dx / 28;
offset = Nx * dx / 3;

% define initial pressure distribution on regular grid as a Gaussian
p0 = exp(-((x - offset) / width).^2);

%% ========================================================================
% RUN SIMULATIONS USING DTT-BASED PSTD METHOD
% =========================================================================

% calculate the time step
dt = CFL * dx / c0;

% normalised plot axes
x_ax = (0:Nx - 1) / (Nx - 1);
t_ax = (0:waterfall_snapshots - 1) / (waterfall_snapshots - 1);

% preallocate matrix to store simulation data
p_four_bc = zeros(bc_num, Nx, waterfall_snapshots);

% loop over boundary conditions
for bc_ind = 1:bc_num
    
    % initialise index for storing waterfall data
    waterfall_ind = 1;

    % define which transforms to use
    switch bc_ind
        case 1
            
            % Neumann-Neumann / WSWS
            T1 = DCT1;
            T2 = DST2;
            
            % set indices for representative sample
            ind1 = 1;
            ind2 = Nx;
            
        case 2
            
            % Neumann-Dirichlet / WSWA
            T1 = DCT3;
            T2 = DST4;
            
            % set indices for representative sample
            ind1 = 1;
            ind2 = Nx - 1;
            
        case 3
            
            % Dirichlet-Neumann / WAWS
            T1 = DST3;
            T2 = DCT4;
            
            % set indices for representative sample
            ind1 = 2;
            ind2 = Nx;

        case 4

            % Dirichlet-Dirichlet / WAWA
            T1 = DST1;
            T2 = DCT2;
            
            % set indices for representative sample
            ind1 = 2;
            ind2 = Nx - 1;

    end
    
    % assign initial condition for p, just selecting representative sample
    p = p0(ind1:ind2);
    
    % run the model backwards for dt/2 to calculate the initial condition
    % for u, to take into account the time staggering between p and u, in
    % this case setting u0 = 0
    u = (dt / 2) / rho0 * gradientDtt1D(p, dx, T1, shift_output, align_output);
        
    % calculate pressure in a loop
    for time_ind = 1:Nt

        % update u
        u = u - dt / rho0 * gradientDtt1D(p, dx, T1, shift_output, align_output);

        % update p
        p = p - dt * rho0 * c0^2 * gradientDtt1D(u, dx, T2, shift_output, align_output);
        
        % save the pressure field
        if ~rem(time_ind, floor(Nt / waterfall_snapshots))
            p_four_bc(bc_ind, ind1:ind2, waterfall_ind) = p;
            waterfall_ind = waterfall_ind + 1;
        end

    end
    
end

%% ========================================================================
% COMBINE SIMULATIONS WITH DIFFERENT BC
% =========================================================================

% form matrix of the different possible boundary conditions, where +1
% corresponds to a positive image source, and -1 to a negative image
% source, thus the four columns correspond to Neumann-Neumann,
% Neumann-Dirichlet, Dirichlet-Neumann, and Dirichlet-Dirichlet
M = [1  1 -1 -1; 
     1  1  1  1; 
     1 -1  1 -1; 
     1 -1 -1  1];

 % open a new figure window
plot_fig = figure;
 
% loop over different reflection coefficients
for ind = 1:2
     
    % set the desired reflection coefficient
    switch ind
        case 1
             
            % non-reflecting boundaries
            R_l = 0;
            R_r = 0;
             
        case 2
             
            % partially reflecting boundary
            R_l = 0;
            R_r = 0.5;
            
    end
     
    % form r vector (this corresponds to the image source amplitudes)
    r = [R_l, 1, R_r, R_r * R_l].';
    
    % solve for the weights for each of the pre-calculated fields
    w = 1/4 * M.' * r %#ok<NOPTS>
     
    % preallocate matrix
    p_waterfall = zeros(Nx, waterfall_snapshots);
     
    % loop over time (waterfall) index
    for waterfall_ind = 1:waterfall_snapshots
        
        % form the pressure field by summing weighted fields
        p_waterfall(:, waterfall_ind) = w.' * p_four_bc(:, :, waterfall_ind);
        
        % plot the fields and the summation
        figure(plot_fig);
        
        subplot(5, 1, 1);
        plot(x_ax, p_four_bc(1, :, waterfall_ind), 'k-');
        set(gca, 'YLim', [-1, 1]);
        title('Neumann-Neumann (WSWS)');
        
        subplot(5, 1, 2);
        plot(x_ax, p_four_bc(2, :, waterfall_ind), 'k-');
        set(gca, 'YLim', [-1, 1]);
        title('Neumann-Dirichlet (WSWA)');
        
        subplot(5, 1, 3);
        plot(x_ax, p_four_bc(3, :, waterfall_ind), 'k-');
        set(gca, 'YLim', [-1, 1]);
        title('Dirichlet-Neumann (WAWS)');
        
        subplot(5, 1, 4);
        plot(x_ax, p_four_bc(4, :, waterfall_ind), 'k-');
        set(gca, 'YLim', [-1, 1]);
        title('Dirichlet-Dirichlet (WAWA)');
        
        subplot(5, 1, 5);
        plot(x_ax, p_waterfall(:, waterfall_ind), 'k-');
        set(gca, 'YLim', [-1, 1]);
        title(['Sum r = [' num2str(r(1)) ', ' num2str(r(2)) ', ' num2str(r(3)) ', ' num2str(r(4)) ']']);
        drawnow;
        
    end
     
    % waterfall plot of evolution of field
    figure;
    waterfall(x_ax, t_ax, p_waterfall.');
    view(10, 70);
    colormap([0, 0, 0]);
    set(gca, 'ZLim', [-1, 1], 'FontSize', 12);
    ylabel('time');
    zlabel('pressure');
    xlabel('position');
    title(['r = [' num2str(r(1)) ', ' num2str(r(2)) ', ' num2str(r(3)) ', ' num2str(r(4)) ']']);
    grid off;
     
end

% close animation figure
close(plot_fig);