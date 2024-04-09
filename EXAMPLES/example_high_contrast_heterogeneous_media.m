% Example code
%
% Simulates scattering of a plane wave by an arrangement of high-contrast 
% penetrable square particles.
%
% This is the application in Section 3c) of [1]. See also Section 6 of [2]
% for more details about scattering by high-contrast heterogeneous media.
%
% References:
%
% [1] Stuart C. Hawkins et al. Metamaterial applications of TMATSOLVER, an
% easy-to-use software for simulating multiple scattering in two dimensions,
% Proceedings of the Royal Society A, 2024.
%
% [2] D Peterseim, B Verfürth. Computational high frequency scattering from 
% high-contrast heterogeneous media, Mathematics of Computation 89 (326), 
% pp. 2649-2674, 2020.
%
%
% Stuart C. Hawkins, Barbara Verfuerth - 6 April 2024

% Copyright 2023, 2024 Stuart C. Hawkins, Luke G. Bennetts, Malte A. Peter and
% Barbara Verfuerth.
% 	
% This file is part of TMATSOLVER.
% 
% TMATSOLVER is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATSOLVER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATSOLVER.  If not, see <http://www.gnu.org/licenses/>.


clear all

%-----------------------------------------------
% set parameters
%-----------------------------------------------

% half the width of the sides of the square scatterers
half_width = 0.25/8;

% select which figure to reproduce
mode = 8;

%-----------------------------------------------
% set pattern matrices
%-----------------------------------------------

switch mode

    case 1

        % Figure 5 (top left) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            0 1 0 0 0
            ];

    case 2

        % Figure 5 (top right) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 1 1 0 0
            ];

    case 3

        % Figure 5 (bottom left) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            ];

    case 4

        % Figure 5 (bottom right) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            ];

    case 5

        % Figure 6 (top left) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 0 0 0 0
            0 0 0 0 0
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 1 1 1 1
            0 0 0 0 0
            0 0 0 0 0
            ];

    case 6

        % Figure 6 (top right) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 0 0 0 1
            0 0 1 0 1
            0 1 1 0 0
            0 0 1 0 1
            0 0 0 0 1
            0 1 1 0 1
            0 1 1 1 1
            0 0 0 1 0
            ];

    case 7

        % Figure 6 (bottom left) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 1 0 1 0
            0 0 1 1 0
            0 0 1 0 0
            0 1 0 0 0
            0 1 1 0 1
            0 1 1 1 0
            0 1 1 0 0
            0 1 0 0 1
            ];

    case 8

        % Figure 6 (bottom right) in [1]

        % scatterers are arranged in a grid and the 1s in pattern mark where the
        % scatterers are placed
        pattern = [
            0 0 0 1 0
            0 1 1 1 0
            0 0 1 0 0
            0 1 1 0 1
            0 1 0 1 0
            0 0 1 1 0
            0 1 0 0 1
            0 0 0 1 1
            ];

    otherwise

        error('mode must be between 1 and 8.')

end

%-----------------------------------------------
% load T-matrix
%-----------------------------------------------

% load the saved T-matrix
T = tmatrix('tmat_square64_k28p00000.mat');

% extract the wavenumber from the T-matrix
kwave = T.kwave;

%-----------------------------------------------
% set dependent parameters
%-----------------------------------------------

% incident plane wave with direction (cos 0,sin 0)
uinc = plane_wave(0,kwave);

% setup the scatterers, using the loaded T-matrix... set the radius 
% for the mask to be sqrt(2)*half_width ie the distance from the origin to
% the corners of the square
type = extern_square(half_width,half_width,kwave,sqrt(2)*half_width,'',T);

% find the nonzeros in the pattern, which indicate where the scatterers
% should be
[j,i] = find(pattern);

% note we need to invert j because indexing into matrices is down the
% columns
j = size(pattern,1)-j;

% convert i and j into (x,y) coordinates of the centers of the scatterers, 
% which are stored as real and imaginary parts of a complex numbebr
pos = (0.5 + i + 0.5i + 1i * j)/8;

%-----------------------------------------------
% setup and solve scattering problem
%-----------------------------------------------

% note that the solver is instantiated with only the incident
% wave... the scatterers are added later
obj = tmatsolver(uinc);

% add the scatterers to the solver
for j=1:length(pos)
    obj.addParticle(particle(type,pos(j)));
end

% solve the scattering problem
obj.solve()

%-----------------------------------------------
% plot the real part of the total field
%-----------------------------------------------

% plot the total field on [0,1.5]x[0,1]
obj.plot([0 1.5 0 1])

% make the plot look nice
axis tight
colorbar
