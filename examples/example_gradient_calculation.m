% DESCRIPTION:
%     This example script creates different periodic sequences of known
%     symmetry, and then numerically calculates the gradients using
%     gradientDtt1D. The numerical gradient is then compared with the known
%     analytical gradient. A test sequence is created for all supported DTT
%     symmetries. The grid shifting functionality of gradientDtt1D can also
%     be tested by changing the value of shift.
%       
% ABOUT:
%     author      - Bradley Treeby
%     date        - 23 April 2013
%     last update - 1 April 2020
%
% Copyright (C) 2013-2020 Bradley Treeby
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

% shift value
%   0: no shift
%   1: shift by + dx/2
%   2: shift by - dx/2
shift = 0;

% set test case
test_cases = 1:8;

% loop through test cases
for index = 1:length(test_cases)
    
    test_case = test_cases(index);

    switch test_case

    % =====================================================================
    % DCT-I WSWS
    % =====================================================================

        case 1
            
            % create function with appropriate symmetry
            Nx              = 10;
            dx              = pi/(Nx-1);
            x_in            = 0:dx:pi;
            f               = cos(x_in);
            
            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end

            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = -sin(x_out);
            
            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = cos(x_full);
            dfdx_full       = -sin(x_full);

            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DCT-I WSWS');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');
            
    % =====================================================================
    % DCT-II HSHS
    % =====================================================================

        case 2

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = pi/(Nx-1);
            x_in            = dx/2:dx:pi-dx/2;
            f               = cos(x_in);

            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            
            
            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = -sin(x_out);

            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = cos(x_full);
            dfdx_full       = -sin(x_full);
            
            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DCT-II HSHS');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DCT-III WSWA
    % =====================================================================

        case 3

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = (pi/2)/(Nx);
            x_in            = 0:dx:pi/2-dx;
            f               = cos(x_in);
            
            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            

            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = -sin(x_out);
            
            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = cos(x_full);
            dfdx_full       = -sin(x_full);

            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DCT-III WSWA');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DCT-IV HSHA
    % =====================================================================

        case 4

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = (pi/2)/(Nx-1);
            x_in            = dx/2:dx:pi/2-dx/2;
            f               = cos(x_in);
            
            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            

            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = -sin(x_out);

            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = cos(x_full);
            dfdx_full       = -sin(x_full);
            
            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DCT-IV HSHA');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DST-I WAWA
    % =====================================================================     

        case 5
            
            % create function with appropriate symmetry
            Nx              = 10;
            dx              = pi/(Nx+1);
            x_in            = dx:dx:pi-dx;
            f               = sin(x_in);

            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            
            
            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = cos(x_out);

            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = sin(x_full);
            dfdx_full       = cos(x_full);
            
            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DST-I WAWA');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DST-II HAHA
    % =====================================================================     

        case 6

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = pi/Nx;
            x_in            = dx/2:dx:pi-dx/2;
            f               = sin(x_in);

            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            
            
            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = cos(x_out);
            
            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = sin(x_full);
            dfdx_full       = cos(x_full);

            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DST-II HAHA');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DST-III WAWS
    % =====================================================================

        case 7

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = (pi/2)/Nx;
            x_in            = dx:dx:pi/2;
            f               = sin(x_in);

            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            
            
            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = cos(x_out);
            
            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = sin(x_full);
            dfdx_full       = cos(x_full);

            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DST-III WAWS');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================
    % DST-IV HAHS
    % =====================================================================

        case 8

            % create function with appropriate symmetry
            Nx              = 10;
            dx              = (pi/2)/(Nx-1);
            x_in            = dx/2:dx:pi/2-dx/2;
            f               = sin(x_in);
            
            switch shift
                case 0
                    x_out = x_in;
                case 1
                    x_out = x_in + dx/2;
                case 2
                    x_out = x_in - dx/2;
            end            

            % compute gradient 
            dfdx_num        = gradientDtt1D(f, dx, test_case, shift);
            dfdx_analytical = cos(x_out);
            
            % plot functions
            x_full          = -dx:dx/5:2*pi - dx;
            f_full          = sin(x_full);
            dfdx_full       = cos(x_full);

            % plot and compare
            figure;
            subplot(2, 1, 1);
            plot(x_full, f_full, 'k-', x_full, dfdx_full, 'b-', x_out, dfdx_num, 'r.', x_in, f, 'k.');
            legend('Input Function', 'Analytical Derivative', 'Numerical Derivative', 'Location', 'Best');
            xlabel('x');
            ylabel('f(x) and f''(x)');
            title('DST-IV HAHS');
            set(gca, 'XLim', [-dx, 2*pi - dx]);
            
            subplot(2, 1, 2);
            plot(x_out, abs(dfdx_num - dfdx_analytical), 'k.-');
            xlabel('x');
            ylabel('error');
            set(gcf,'color','w');

    % =====================================================================     

    end
    
end