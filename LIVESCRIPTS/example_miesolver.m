% Example code
%
% Simulates scattering by equally spaced circular-cylinders arranged on a
% circle.
%
% Compares the far field with that computed by MieSolver [1].
%
% References:
%
% [1] Stuart C. Hawkins, ACM Trans. Math. Softw. Vol 46 Article 19 (2020)
%
% Stuart C. Hawkins - 5 April 2024

% Copyright 2023, 2024 Stuart C. Hawkins.
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

%-----------------------------------------
% check for MieSolver
%-----------------------------------------

% MieSolver is a package for scattering by circular cylinders.
% It is available at
%
% http://www.miesolver.org
%
% Please see Stuart C. Hawkins, ACM Trans. Math. Softw. Vol 46 Article 19 (2020)
% for full details.

% if miesolver.m does not exist as a file then probably MieSolver
% is not installed.
if ~exist('miesolver.m','file')
    error('Please download and install MieSolver from www.miesolver.org')
end

%-----------------------------------------
% set parameters
%-----------------------------------------

% wavenumber
kwave = 2*pi;

% scatterers lie on a circle with radius R
R = 10;

% number of scatterers
N = 10;

%-----------------------------------------
% set dependent parameters
%-----------------------------------------

% work out positions of scatterers evenly spaces on the circle of radius R
pos = R*exp(1i*2*pi*(0:N-1)/N);

% set radii of scatterers (all have unit radius)
rad = ones(size(pos));

% set incident wave
inc = plane_wave(0,kwave);

% setup vector of observation directions for the far field
tp = linspace(0,2*pi,1000);
z = exp(1i*tp);

%-----------------------------------------
% Mie version
%-----------------------------------------

% setup MieSolver object to perform the calculation
pr = MieSolver(inc);

% add the scatterers (as sound soft circular cylinders)
for j=1:length(pos)
    pr.addScatterer(scatterer(pos(j),rad(j),'SOFT'));
end

% solve the scattering problem
pr.solve();

% get the far field
f = pr.getFarfield(z);

%-----------------------------------------
% TMATSOLVER version
%-----------------------------------------

% get the order for the wavefunction expansions
n = suggestedorder(kwave,rad(1));

% set the scatter types (only need one... a sound soft circle)
types = [sound_soft_disk(kwave,rad(1))];

% setup the TMATSOLVER object to perform the calculation
tm = tmatsolver(inc);

% adcd the scatterers - all are type 1
for j=1:length(pos)
    tm.addParticle(particle(types(1),pos(j)));
end

% solve the scattering problem
tm.solve()

% get the far field
g = tm.getFarfield(z);

%-----------------------------------------
% compute error
%-----------------------------------------

% compute and print the error
fprintf('Error in farfield between TMATSOLVER and MieSolver: %0.2e\n',max(abs(f-g)));

%-----------------------------------------
% figure to show the configuration
%-----------------------------------------

figure(1)

% vector of angles that we need to plot the circle of radius R
phi = linspace(0,2*pi,100);

% plot the circle of radius R on which the scatterers sit
plot(R*cos(phi),R*sin(phi),'-','color',0.8*[1 1 1]);

% add the scatterers
hold on
tm.schematic

% make the figure look nice
title('Scattering configuration')
axis equal
hold off

%-----------------------------------------
% plot the far fields
%-----------------------------------------

figure(2)

% plot the far field in dB
plot(tp,10*log10(2*pi*abs(f).^2),'r-',...
    tp,10*log10(2*pi*abs(g).^2),'b-')

% add a legend
legend('MieSolver','TMATSOLVER')

% add a title
title('scattering cross section (dB)')

