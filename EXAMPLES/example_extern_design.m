% Demonstration of multiple scattering using a T-matrix loaded from outside
%
% Stuart C. Hawkins - 21 June 2023

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

%-----------------------------------------------
% set parameters
%-----------------------------------------------

% separation of centres of the scatterers
d = 4;

% filename of T-matrix
fname = 'saved_tmatrix.mat';

% number of scatterers
N = 40;

% radius of scatterers
% Note: we don't always know this... that's not too important, it is just
% used to mask out the scatterers when we do a plot, it isn't used for any
% calculations. We can just put a dummy value
rad = sqrt(2);

%-----------------------------------------------
% load the T-matrix
%-----------------------------------------------

% load the T-matrix and construct tmatrix object
T = tmatrix(fname);

%-----------------------------------------------
% get parameters from the T-matrix
%-----------------------------------------------

% wavenumber
kwave = T.kwave;

% truncation parameter
n = T.order;

%-----------------------------------------------
% setup experiment
%-----------------------------------------------

% positions of scatterers
pos = d*(1:N);

% incident field (point source)
uinc = point_source(-d,kwave);

% setup array of types containing the scatterer
types = [extern_design(kwave,rad,'',T)];

%-----------------------------------------------
% setup and solve scattering problem
%-----------------------------------------------

% note that the solver is instantiated with only the incident 
% wave... the scatterers are added later
obj = tmatsolver(uinc);

% add the scatterers to the solver
for j=1:length(pos)
    obj.addParticle(particle(types(1),pos(j)));
end

% solve the scattering problem
obj.solve()

% setup a rectange specifying the four corners of the box in which we want
% to plot the field
R = [pos(1)-10 pos(end)+10 -4*d 4*d];

% plot the total field
obj.plot(R)

% % make 
% axis(R)
