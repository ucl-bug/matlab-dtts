function deriv = gradientDtt1D(f, dx, dtt_type, shift, align_output)
%GRADIENTDTT1D Calculate gradient using discrete trigonometric transforms.
%
% DESCRIPTION:
%     gradientDtt1D computes the spectral gradient of a 1D input function
%     using discrete trigonometric transforms (DTTs). The DTTs used to
%     transform the input to and from the frequency domain are chosen based
%     on the symmetry of the input f, defined using dtt_type. The gradient
%     values can also be returned on a staggered grid. 
%
%     Note, the group I DTTs (DCT-I, DCT-II, DST-I, DST-II), have different
%     length k-space vectors due to the implied symmetry of the input
%     sequence (e.g., the Nyquist and/or DC components are zero, and do not
%     need to be defined). This means that in some cases the basis
%     functions weights are trimmed or appended after the forward
%     transform. To ensure the output of gradientDTT1D is the same length
%     as the input f, additional values can be appended or removed from the
%     calculated gradient after the inverse transform is calculated by
%     setting the optional input pad to true (the default). These values
%     are known from the symmetry of the output sequence.
%
%     For additional details on gradient calculation using DTTs, see [1].
%
%     [1] E. Wise, J. Jaros, B. Cox, and B. Treeby, "Pseudospectral
%     time-domain (PSTD) methods for the wave equation: Realising boundary
%     conditions with discrete sine and cosine transforms", 2020.
%
% INPUTS:
%     f            - Vector to find the gradient of.
%     dx           - Grid point spacing.
%     dtt_type     - Type of discrete trigonometric transform. This should
%                    correspond to the assumed input symmetry of the input
%                    function, where:
%
%                        1: DCT-I    WSWS
%                        2: DCT-II   HSHS
%                        3: DCT-III  WSWA
%                        4: DCT-IV   HSHA
%                        5: DST-I    WAWA
%                        6: DST-II   HAHA
%                        7: DST-III  WAWS
%                        8: DST-IV   HAHS
%
% OPTIONAL INPUTS:
%     shift        - Integer controlling whether derivative is shifted to a
%                    staggered grid (default = 0), where 
%
%                        0: no shift
%                        1: shift by + dx/2
%                        2: shift by - dx/2
%
%     align_output - Boolean controlling whether the returned values are
%                    padded and trimmed based on the implied symmetry so
%                    the output is the same length as the input (default =
%                    true). Note, if align_output is false, then the
%                    gradient calculated for shift = 1 and shift = 2 will
%                    be the same. 
%
% OUTPUTS:
%     dfdx         - Gradient of the input function.
%       
% ABOUT:
%     author       - Bradley Treeby
%     date         - 23 April 2013
%     last update  - 20 April 2020
%
% Copyright (C) 2013-2020 Bradley Treeby
%
% See also dtt1D

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

% define literals for DTT types used by the function dtt1D 
DCT1 = 1;    % WSWS
DCT2 = 2;    % HSHS
DCT3 = 3;    % WSWA
DCT4 = 4;    % HSHA
DST1 = 5;    % WAWA
DST2 = 6;    % HAHA
DST3 = 7;    % WAWS
DST4 = 8;    % HAHS

% check for shift input
if (nargin < 4) || isempty(shift)
    shift = 0;
end

% check for pad input
if (nargin < 5) || isempty(align_output)
    align_output = true;
end

% check inputs
validateattributes(f,            {'numeric'}, {'real', 'vector', 'nonsparse'},                      'gradientDTT1D', 'f',            1);
validateattributes(dx,           {'numeric'}, {'real', 'scalar', 'nonsparse'},                      'gradientDTT1D', 'dx',           2);
validateattributes(dtt_type,     {'numeric'}, {'integer', 'scalar', 'nonsparse', '>=', 1, '<=', 8}, 'gradientDTT1D', 'dtt_type',     3);
validateattributes(shift,        {'numeric'}, {'integer', 'scalar', 'nonsparse', '>=', 0, '<=', 2}, 'gradientDTT1D', 'shift',        4);
validateattributes(align_output, {'logical'}, {'scalar'},                                           'gradientDTT1D', 'align_output', 5);

% reshape to be a row vector
f  = reshape(f, 1, []);

% extract length of input function
Nx = length(f);

% compute the implied period of the input function
switch dtt_type
    case 1
        M = 2 .* (Nx - 1);
    case 5
        M = 2 .* (Nx + 1);
    otherwise
        M = 2 .* Nx;
end

% calculate the wavenumbers
switch dtt_type
    case 1

        % WSWS / DCT-I
        n = (0:1:M/2);
        kx = 2 .* pi .* n ./ (M .* dx);

    case 2

        % HSHS / DCT-II
        n = (0:1:(M/2 - 1));
        kx = 2 .* pi .* n ./ (M .* dx);

    case 5

        % WAWA / DST-I
        n = (1:1:(M/2 - 1));
        kx = 2 .* pi .* n ./ (M .* dx);

    case 6

        % HAHA / DST-II
        n = (1:1:M/2);
        kx = 2 .* pi .* n ./ (M .* dx);

    case {3, 4, 7, 8}

        % WSWA / DCT-III
        % HSHA / DCT-IV
        % WAWS / DST-III
        % HAHS / DST-IV
        n = (0:1:(M/2 - 1));
        kx = 2 .* pi .* (n + 0.5) ./ (M .* dx);

end

% compute forward transform and multiply by the wavenumbers
if dtt_type < 5
    
    % discrete cosine transform
    deriv = -kx .* dtt1D(f, dtt_type);
    
else
    
    % discrete sin transform
    deriv =  kx .* dtt1D(f, dtt_type);
    
end

% select appropriate function range for inverse transform
% (this is only required for the whole-wavenumber DTTs as the function
% range is the same for the half-wavenumber DTTs) 
switch dtt_type
    case 1
            
        if shift
            
            % WSWS -> HAHA, remove left endpoint
            deriv = deriv(2:end);
            
        else
        
            % WSWS -> WAWA, remove both endpoints
            deriv = deriv(2:end - 1);
            
        end
        
    case 2
        
        if shift
            
            % HSHS -> WAWA, remove left endpoint
            deriv = deriv(2:end);
            
        else
        
            % HSHS -> HAHA, remove left endpoint and append a zero
            deriv = [deriv(2:end), 0];
            
        end
        
    case 5
        
        if shift
            
            % WAWA -> HSHS, prepend zero
            deriv = [0, deriv];
            
        else
        
            % WAWA -> WSWS, prepend and append zeros
            deriv = [0, deriv, 0];
            
        end
        
    case 6
        
        if shift
            
            % HAHA -> WSWS, prepend zero
            deriv = [0, deriv];
            
        else
        
            % HAHA -> HSHS, prepend zero, remove right endpoint
            deriv = [0, deriv(1:end - 1)];
            
        end
        
end

% compute inverse transform and normalise by the implied period of the
% input function
switch dtt_type
    case 1
        
        if shift
            
            % WSWS -> HAHA, S2^-1 = S3
            deriv = dtt1D(deriv, DST3) ./ M;
            
        else
        
            % WSWS -> WAWA, S1^-1 = S1
            deriv = dtt1D(deriv, DST1) ./ M;
            
        end
        
    case 2
        
        if shift
            
            % HSHS -> WAWA, S1^-1 = S1
            deriv = dtt1D(deriv, DST1) ./ M;
            
        else
        
            % HSHS -> HAHA, S2^-1 = S3
            deriv = dtt1D(deriv, DST3) ./ M;
            
        end
        
    case 3
        
        if shift
            
            % WSWA -> HAHS, S4^-1 = S4
            deriv = dtt1D(deriv, DST4) ./ M;
            
        else
        
            % WSWA -> WAWS, S3^-1 = S2
            deriv = dtt1D(deriv, DST2) ./ M;
            
        end
        
    case 4
        
        if shift
            
            % HSHA -> WAWS, S3^-1 = S2
            deriv = dtt1D(deriv, DST2) ./ M;   
            
        else
            
            % HSHA -> HAHS, S4^-1 = S4
            deriv = dtt1D(deriv, DST4) ./ M;
            
        end
        
        
    case 5
        
        if shift
            
            % WAWA -> HSHS, C2^-1 = C3
            deriv = dtt1D(deriv, DCT3) ./ M;
            
        else
        
            % WAWA -> WSWS, C1^-1 = C1
            deriv = dtt1D(deriv, DCT1) ./ M;
            
        end
        
    case 6
        
        if shift
            
            % HAHA -> WSWS, C1^-1 = C1
            deriv = dtt1D(deriv, DCT1) ./ M;  
            
        else
            
            % HAHA -> HSHS, C2^-1 = C3
            deriv = dtt1D(deriv, DCT3) ./ M;            
            
        end
        
    case 7
        
        if shift
            
            % WAWS -> HSHA, C4^-1 = C4
            deriv = dtt1D(deriv, DCT4) ./ M;
            
        else

            % WAWS -> WSWA, C3^-1 = C2
            deriv = dtt1D(deriv, DCT2) ./ M;
            
        end
        
    case 8
        
        if shift

            % HAHS ->WSWA, C3^-1 = C2
            deriv = dtt1D(deriv, DCT2) ./ M;            
            
        else
            
            % HAHS -> HSHA, C4^-1 = C4
            deriv = dtt1D(deriv, DCT4) ./ M;     
            
        end
        
end

% add back in the implied values so the output is the same length as the
% input
if align_output
    switch dtt_type
        case 1

            switch shift
                case 0

                    % WSWS -> WAWA, add both endpoints
                    deriv = [0, deriv, 0];

                case 1

                    % WSWS -> HAHA, shift right, mirror right endpoint
                    deriv = [deriv, -deriv(end)];                

                case 2

                    % WSWS -> HAHA, shift left, mirror left endpoint
                    deriv = [-deriv(1), deriv];                

            end

        case 2

            switch shift
                case 0

                    % HSHS -> HAHA, no change

                case 1

                    % HSHS -> WAWA, shift right, append zero
                    deriv = [deriv, 0];

                case 2

                    % HSHS -> WAWA, shift left, prepend zero
                    deriv = [0, deriv];

            end

        case 3

            switch shift
                case 0

                    % WSWA -> WAWS, prepend zero, remove right endpoint
                    deriv = [0, deriv(1:end - 1)];                

                case 1

                    % WSWA -> HAHS, shift right, no change

                case 2

                    % WSWA -> HAHS, shift left, mirror left endpoint, remove
                    % right endpoint
                    deriv = [-deriv(1), deriv(1:end - 1)];

            end

        case 4

            switch shift
                case 0

                    % HSHA -> HAHS, no change

                case 1

                    % HSHA -> WAWS, shift right, no change

                case 2

                    % HSHA -> WAWS, shift left, prepend zero, remove right
                    % endpoint
                    deriv = [0, deriv(1:end - 1)];

            end

        case 5

            switch shift
                case 0

                    % WAWA -> WSWS, remove both endpoints
                    deriv = deriv(2:end - 1);

                case 1

                    % WAWA -> HSHS, shift right, remove left endpoint
                    deriv = deriv(2:end);

                case 2

                    % WAWA -> HSHS, shift left, remove right endpoint
                    deriv = deriv(1:end - 1); 

            end

        case 6

            switch shift
                case 0

                    % HAHA -> HSHS, no change

                case 1

                    % HAHA -> WSWS, shift right, remove left endpoint
                    deriv = deriv(2:end);

                case 2

                    % HAHA -> WSWS, shift left, remove right endpoint
                    deriv = deriv(1:end - 1); 

            end

        case 7

            switch shift
                case 0

                    % WAWS -> WSWA, remove left endpoint, append zero
                    deriv = [deriv(2:end), 0];

                case 1

                    % WAWS -> HSHA, shift right, remove left endpoint, mirror
                    % right endpoint 
                    deriv = [deriv(2:end), -deriv(end)];

                case 2

                    % WAWS -> HSHA, shift left, no change

            end

        case 8

            switch shift
                case 0

                    % HAHS -> HSHA, no change

                case 1

                    % HAHS -> WSWA, shift right, remove left endpoint, append
                    % zero 
                    deriv = [deriv(2:end), 0];

                case 2

                    % HAHS -> WSWA, shift left, no change

            end

    end
end