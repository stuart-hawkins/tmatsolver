% Demonstration of multiple scattering using TMATSOLVER
%
% Stuart C. Hawkins - 10 January 2023

% Copyright 2023 Stuart C. Hawkins.
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
close all

%-----------------------------------------------
% set parameters
%-----------------------------------------------

% incident wavelength
lambda = 2;

% radius of scatterers
rho = 1;

% number of scatterers
N = 8;

% radius of circle on which scatterers placed
R = 10;

%-----------------------------------------------
% set dependent parameters
%-----------------------------------------------

% wavenumber
kwave = 2*pi/lambda;

% incident plane wave with direction (cos 0,sin 0)
uinc = plane_wave(0,kwave);

% setup types of scatterers 
types = [sound_soft_disk(kwave,rho),...
    sound_hard_disk(kwave,rho)];

% place the scatterers equally spaced on circle
pos = R*exp(1i*2*pi*(0:N-1)/N);

%-----------------------------------------------
% setup and solve scattering problem
%-----------------------------------------------

% note that the solver is instantiated with only the incident 
% wave... the scatterers are added later
obj = tmatsolver(uinc);

% add the scatterers to the solver
for j=1:length(pos)
    % use alternating type ie sound soft/sound hard
    obj.addParticle(particle(types(rem(j+1,2)+1),pos(j)));
end

% solve the scattering problem
obj.solve()

%-----------------------------------------------
% plot movie
%-----------------------------------------------

% half-width of square on which to do the plot
L = 2*R;

% generate the movie
m = obj.movie([-L L -L L]);

% export the movie
m.export('movie.gif','gif')
fprintf('Movie written to movie.gif\n');
