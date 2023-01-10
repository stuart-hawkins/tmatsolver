% Compute T-matrix of a circular scatterer
%
%  T = mie_tmatrix(k,r,n,'SOFT') computes the T-matrix for a sound-soft
%  circular scatterer with radius r. Here k is the wavenumber and n is the
%  order of the T-matrix.
%
%  T = mie_tmatrix(k,r,n,'HARD') computes the T-matrix for a sound-hard
%  circular scatterer with radius r. Here k is the wavenumber and n is the
%  order of the T-matrix.
%
%  T is of type tmatrix (see the TMATROM package http://www.romapp.org). 
%
% See also: tmatrix.
%
% Stuart C. Hawkins - 7 January 2023

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


function T = mie_tmatrix(kwave,radius,n,BC,x0)

% set default for scatterer centre
if nargin<5
    x0 = 0;
end

% precompute Hankel and Bessel function values
l = -n:n;
h = besselh(0:n,kwave*radius);
j = besselj(0:n,kwave*radius);

% ...and their derivatives
hd = besselhd(0:n,kwave*radius);
jd = besseljd(0:n,kwave*radius);

% H and J contain values of Hankel and Bessel functions for orders
% -l:l
H = h(abs(l)+1);
J = j(abs(l)+1);

% ...and similar for derivatives
Hd = hd(abs(l)+1);
Jd = jd(abs(l)+1);

% compute the T-matrix entries... the T-matrices are diagonal because
% the scatterers are circular
switch upper(BC)

    case 'SOFT'

        M = diag(-J./H);

    case 'HARD'

        M = diag(-Jd./Hd);

    otherwise

        error('Boundary condition %s not recognised/supported',BC)

end

% create a T-matrix object representing the T-matrix
T = tmatrix(n,kwave,M,x0);

