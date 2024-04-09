% Class to represent 'types' of scatterer in tmatsolver.
%
%   The main idea is that this class associates details about the scatterer
%   with its T-matrix, so that eg we can plot the scatterer.
%
% Stuart C. Hawkins - 7 January 2023

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


classdef design < handle & matlab.mixin.Heterogeneous

    % Note: matlab.mixin.Heterogeneous means we can have arrays of this
    % class and its children

    properties
        name
        id
        tmat
    end

    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = design(tmat,name)            

            % set default for name
            if nargin<2
                name = '';
            end

            % store the T-matrix
            self.tmat = tmat;

            % store the name
            self.name = name;

            % generate and store an id number
            self.id = design.generate_id();
   
        end

        %-----------------------------------------------
        % T-matrix error measure
        %-----------------------------------------------

        function err = error(self)

            % use T-matrix error method
            err = self.tmat.error();

        end

        %-----------------------------------------------
        % default plot function... this can be overloaded
        % in child classes
        %-----------------------------------------------

        function varargout = plot(self)

            % do nothing

            % return an empy handle if needed
            if nargout>0
                varargout{1} = [];
            end

        end

    end

    %=================================================================
    % static methods
    %=================================================================

    methods(Static)

        %-----------------------------------------------
        % generate an id number
        %-----------------------------------------------

        function val = generate_id()

            % static variable is implemented in Matlab using a persistent
            % variable
            persistent design_id_counter

            % initialise the counter if needed
            if isempty(design_id_counter)
                design_id_counter = 0;
            end

            % increment the counter because we are generating a new id
            design_id_counter = design_id_counter + 1;

            % return the new id
            val = design_id_counter;

        end

    end

end