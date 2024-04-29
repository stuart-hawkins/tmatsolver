% Design representing a square scatterer
%
%        scat = extern_square(half_width,half_height,k,r,'',M) sets up an
%        object representing a rectangular scatterer with dimensions
%        2*half_width x 2*half_height. The T-matrix of the scatterer is M. 
%        The wavenumber is k. r is the radius of the circular mask used for
%        plotting the scattered and total fields.
%
%   scat.plot() plots the scatterer.
%
%   h = scat.plot() plots the scatterer and returns the handle h of the
%   graphics object.
%
% See also: design, extern_design, tmatsolver.
%
% Stuart C. Hawkins, Luke Bennetts, Malter Peter - 5 April 2024

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


classdef extern_square < extern_design
    
    properties
        half_width
        half_height
    end
    
    methods
        
        %-----------------------------------------------
        % constructor
        %-----------------------------------------------
        
        function self = extern_square(half_width,half_height,varargin)
           
            self = self@extern_design(varargin{:});
            
            % store the square parameters
            self.half_height = half_height;
            self.half_width = half_width;
            
        end
        
        %-----------------------------------------------
        % plot function
        %-----------------------------------------------

        function varargout = plot(self)

            % extract the dimensions of the rectangle
            w = self.half_width;
            h = self.half_height;
            
            % plot the scatterer
            h = fill([-w w w -w -w],[-h -h h h -h],'k-');

            % return graphics handle if needed
            if nargout>0
                varargout{1} = h;
            end

        end

    end
    
end