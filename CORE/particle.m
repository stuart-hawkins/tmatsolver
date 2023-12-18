% Class to represent a particle
%
%  In tmsolver particles associate a 'type' of particle with a position.
%
%  p = particle(t,x) represents a particle of type t located at x. t must
%  be of class design.
%
%  p = particle(t,x,phi) represents a particle of type t located at x,
%  rotated by an angle phi in the positive direction. t must be of class
%  design.
%
%  p.plot() visualises the particle. This uses the plot method of the
%  underlying type.
%
%  T = p.tmatrix() returns the T-matrix associated with the particle
%
% See also: tmsolver, design.
%
% Stuart C. Hawkins - 8 January 2023

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


classdef particle < handle

    properties
        type
        pos
        rot
    end

    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = particle(type,pos,rot)

            % set defaults for type, pos and rot in case no arguments are 
            % given... this allows us to create arrays of particle using 
            % the particle.empty construction
            if nargin==0
                type = [];
                pos = NaN;
                rot = NaN;
                return
            end

            % set default for rot
            if nargin<3
                % default is no rotation
                rot = 0;
            end

            % check that type is of the correct class
            if ~isa(type,'design')
                error('type must be of class design')
            end

            % store type of particle
            self.type = type;

            % store position of particle
            self.pos = pos;

            % store orientation of particle
            self.rot = rot;

        end

        %-----------------------------------------------
        % get T-matrix
        %-----------------------------------------------

        function T = tmatrix(self)

            % return the T-matrix of the underlying type
            T = self.type.tmat;

        end

        %-----------------------------------------------
        % plot
        %-----------------------------------------------

        function varargout = plot(self)

            % use the type's plot method to plot the scatterer at the
            % origin
            h = self.type.plot();

            % use hgttransform to rotate the shape and move it to the right 
            % place...
            hgt = hgtransform;
            for k=1:length(h)
                set(h(k),'parent',hgt);
            end
            M1 = makehgtform('zrotate',self.rot);
            M2 = makehgtform('translate',[real(self.pos);imag(self.pos);0]);
            set(hgt,'Matrix',M2*M1)      

            % return graphics handle if needed
            if nargout>0
                varargout{1} = h;
            end

        end

    end

end
