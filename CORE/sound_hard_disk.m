% Class representing sound hard circular scatterers
%
%   scat = sound_hard_disk(k,r) sets up object for a circle with radius r.
%   The wavenumber is k.
%
%   scat.plot() plots the scatterer.
%
%   h = scat.plot() plots the scatterer and returns the handle h of the
%   graphics object.
%
% See also: design, sound_soft_disk, tmatsolver.
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


classdef sound_hard_disk < design
    
    properties
        radius
    end

    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = sound_hard_disk(kwave,radius,name)

            % set default name
            if nargin<3
                name = 'sound hard disk';
            end

            % get order for T-matrix using Wiscombe's formula
            n = suggestedorder(kwave,radius);

            % compute the T-matrix
            T = mie_tmatrix(kwave,radius,n,'HARD');

            % call the parent constructor to setup the scatterer
            self = self@design(T,name);

            % radius is stored only in case we want to plot the scatterer
            self.radius = radius;

        end

        %-----------------------------------------------
        % plot function
        %-----------------------------------------------

        function varargout = plot(self)

            % get angles in [0,2*pi] to plot the circle
            theta = linspace(0,2*pi,200);

            % plot the circle
            h = plot(self.radius*cos(theta),self.radius*sin(theta),'k-');

            % return graphics handle if needed
            if nargout>0
                varargout{1} = h;
            end

        end

    end

end
            

