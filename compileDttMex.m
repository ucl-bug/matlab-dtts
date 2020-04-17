%COMPILEDTTMEX Compile mex-functions for dtt1D, dtt2D, and dtt3D.
%
% DESCRIPTION:
%     compileDttMex compiles the mex functions for dtt1D, dtt2D, and dtt3D.
%
%     On Windows, the compilation uses a pre-compiled version of FFTW
%     included in the repository. If re-compiling, .lib files can be
%     generated from .dll files using a Visual Studio command prompt
%
%         lib /machine:x64 /def:libfftw3-3.def
%     
%     On Linux or macOS, fftw can be downloaded from http://www.fftw.org/.
%     Compile using the following options:
%
%         sudo ./configure --enable-avx CFLAGS="-m64 -fPIC"
%         sudo make
%         sudo make install
% 
%     Note, use --enable-sse2 if avx instructions aren't supported on
%     your processor.
%
%     It is assumed that a suitable C++ compiler is installed and selected
%     by calling mex -setup. See: https://www.mathworks.com/support/...
%     requirements/supported-compilers.html for more details. On linux,
%     specify the version of the compiler and FFTW to use in the code
%     below. On macOS, compiling mex files requires Xcode, which can be
%     downloaded from the app store.   
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 28 September 2017
%     last update   - 19 March 2020
%
% Copyright (C) 2017-2020 Bradley Treeby
%
% See also dtt1D, dtt2D, dtt3D

% check for windows, mac, or linux
if ispc
    
    % use default compiler and link to pre-compiled FFTW library
    mex -L"./" -llibfftw3-3 dtt1D.cpp
    mex -L"./" -llibfftw3-3 dtt2D.cpp
    mex -L"./" -llibfftw3-3 dtt3D.cpp
    
elseif ismac
    
    % use default compiler and link to FFTW installed on local machine
    mex -I/usr/local/include/ -L/usr/local/lib -lfftw3 -lm dtt1D.cpp
    mex -I/usr/local/include/ -L/usr/local/lib -lfftw3 -lm dtt2D.cpp
    mex -I/usr/local/include/ -L/usr/local/lib -lfftw3 -lm dtt3D.cpp

else
    
    % make user update paths based on examples below
    error('Please update paths to GCC and FFTW before calling compileDttMex.');
    
%     % specify gcc compiler and dynamically link to FFTW (change the path
%     % locations in the example below)
%     mex -v GCC='/usr/bin/gcc-7' -I/usr/local/include/ -L/mnt/Apps/software/FFTW/3.3.8-OpenMPI-4.0.1-GCC-7.3.0-2.30/lib -lfftw3 -lm dtt1D.cpp
%     mex -v GCC='/usr/bin/gcc-7' -I/usr/local/include/ -L/mnt/Apps/software/FFTW/3.3.8-OpenMPI-4.0.1-GCC-7.3.0-2.30/lib -lfftw3 -lm dtt2D.cpp
%     mex -v GCC='/usr/bin/gcc-7' -I/usr/local/include/ -L/mnt/Apps/software/FFTW/3.3.8-OpenMPI-4.0.1-GCC-7.3.0-2.30/lib -lfftw3 -lm dtt3D.cpp

end