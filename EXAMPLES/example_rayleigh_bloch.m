% Example code
%
% Simulates scattering by a line array of square cylinders, dominated by a
% Rayleigh-Bloch wave.
%
% This is the application in Section 3d) of [1].
%
% References:
%
% [1] Stuart C. Hawkins et al. Metamaterial applications of TMATSOLVER, an
% easy-to-use software for simulating multiple scattering in two dimensions,
% Proceedings of the Royal Society A, 2024.
%
% Stuart C. Hawkins, Luke G. Bennetts, Malter A. Peter - 6 April 2024

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
% set the main parameters
%-----------------------------------------------

% half the width of the sides of the square scatterers
half_width=0.75;

% number of scatterers
N = 21;

% distance between the scatterers' centres
R = 4;

% set positions of the scatterers
pos = R*[1:N];

% set the box in which we plot the total field
mybox = [pos(1)-10 pos(end)+10 -16*half_width 16*half_width];

%-----------------------------------------------
% load T-matrix
%-----------------------------------------------

% load the saved T-matrix
T = tmatrix('tmatrix66.mat');

% extract the wavenumber from the T-matrix
kwave = T.kwave;

%-----------------------------------------------
% set dependent parameters
%-----------------------------------------------

% incident point source located at the origin
uinc = point_source(0,kwave);

% setup the scatterers, using the loaded T-matrix... set the radius 
% a little smaller so the mask doesn't wipe out the near field
type = extern_square(half_width,half_width,kwave,0.5*half_width,'',T);

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
% plot the total field
%-----------------------------------------------

% evaluate the total field at the midpoints of the scatterers
u = obj.getTotalField(0.5*(pos(1:end-1)+pos(2:end)));

% get the maximum value... we can use this to fix the color axis
mxval = max(abs(u));

% plot the total field
h = obj.plot(mybox);

%-----------------------------------------------
% rescale the x and y axis data to dimensionless
% units x/R and y/R
%-----------------------------------------------

% loop through all handles in h
for j=1:length(h)

    % setup a transform that scales by 1/R
    H = hgtransform;
    M = makehgtform('scale',[1/R 1/R 1]);

    % the scatterers are plotted using translations implemented with 
    % hgtransforms so we can't just scale by setting the transform matrix 
    % to M, that will undo the translation. Instead we get the original
    % transform matrix, and then multiply it by M
    p = get(h(j),'parent');
    if isa(p,'matlab.graphics.primitive.Transform')
        N = get(p,'Matrix');
    else
        % if p isn't a transform then there is no hgtransform already
        % applied so we can just take N to be the identity
        N = eye(4,4);
    end

    % set the new transform matrix
    set(h(j),'parent',H);
    set(H,'Matrix',M*N);
end
    
%-----------------------------------------------
% make the plot look nice
%-----------------------------------------------

% change the axis box to rescaled units
axis(mybox/R)

% make the plot look nice
caxis(mxval*[-1,1])
colorbar
xlabel('x/R')
ylabel('y/R')


