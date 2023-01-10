% Class for solving multiple scattering problems
%
%  obj = tmatsolver(uinc) creates an instance of the class for solving
%  scattering problems with incident wave uinc. Here uinc is of class
%  incident eg of type plane_wave.
%
%  obj.addParticle(p) adds a particle to the scattering problem. Here p is
%  of class particle.
%
%  obj.solve() solves the scattering problem.
%
%  f = obj.getFarfield(x) evaluates the far field f at directions given 
%  in the vector x. The directions are specified by complex numbers with 
%  unit modulus.
%
%  f = obj.getScatteredField(x) evaluates the scattered field at points
%  given in the array x. The points are specified by complex numbers.
%
%  f = obj.getTotalField(x) evaluates the total field at points
%  given in the array x. The points are specified by complex numbers.
%
%  obj.schematic() visualises the scatterers.
%
%  obj.plot(ax) plots the real part of the total field on the domain
%  [ax(1) ax(2)] x [ax(3) ax(4)].
%
%  obj.plot(ax,phase) plots real(phase*u) where u is the total field.
%
%  m = obj.movie(ax) generate a movie of the real part of the total field 
%  on the domain [ax(1) ax(2)] x [ax(3) ax(4)]. The movie is stored in an 
%  object of type mymovie.
%
% See also: particle, plane_wave, mymovie.
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


classdef tmatsolver < handle

    properties
        incident_field
        particles
        soln
    end

    methods

        %----------------------------------
        % constructor
        %----------------------------------

        function self = tmatsolver(incident_field,varargin)

            % store incident field
            self.incident_field = incident_field;

            % other arguments are particles... add them one by one

            % initialise array to empty... this is important because it
            % allows us to add particles later. If we don't initialise the
            % array to have class particle then it will default to real and
            % we won't be allowed to add particles to it later
            self.particles = particle.empty;

            % loop through varargin
            for j=1:length(varargin)

                % check type is correct
                if ~isa(varargin{j},'particle')
                    error('Number %d argument must be of class particle',j+1);
                end

                % check wavenumber of T-matrix matches incident wave
                if abs(self.incident_field.kwave-varargin{j}.tmatrix.kwave) > 1e-12
                    error('T-matrix wavenumber for number %d argument doesn''t match incident field',j+1);
                end

                % add particle to list
                self.particles(j) = varargin{j};

            end

        end

        %----------------------------------
        % solve function
        %----------------------------------

        function solve(self)

            % check incident wave is set
            if isempty(self.incident_field)
                error('First set incident field')
            end

            %- - - - - - - - - - - - - - - - -
            % work out order of each T-matrix
            % and size of matrix etc
            %- - - - - - - - - - - - - - - - -

            % initialise vector holding order of each T-matrix
            orders = zeros(length(self.particles),1);

            % pull out orders
            for j=1:length(self.particles)
                orders(j) = self.particles(j).tmatrix.order;
            end

            % work out degrees of freedom of matrix
            dof = 2*sum(orders)+1;

            % initialise vector holding centres of the particles
            pos = zeros(length(self.particles),1);

            % pull out centres
            for j=1:length(self.particles)
                pos(j) = self.particles(j).pos;
            end

            %- - - - - - - - - - - - - - - - -
            % setup RHS
            %- - - - - - - - - - - - - - - - -

            % loop through the particles
            for j=1:length(self.particles)

                % get wavefunction expansions of the incident field
                a{j} = regularwavefunctionexpansion(orders(j),...
                    pos(j),self.incident_field);

                % temporarily set the origin of the T-matrix so it is
                % compatible with the wavefunction expansion
                self.particles(j).tmatrix.setOrigin(pos(j));

                % apply the T-matrix to the incident field
                b{j} = self.particles(j).tmatrix * a{j};

            end

            % pack the coefficients of b into a vector
            rhs = pack(b);

            %- - - - - - - - - - - - - - - - -
            % solve linear system
            %- - - - - - - - - - - - - - - - -

            % set number of GMRES iterations
            nitns = min(100,floor(dof/2));

            % solve linear system using GMRES
            [x,flag,relres,iter,reshist] = gmres(@matrix_product,rhs,nitns,1e-8,1);

            % check GMRES worked okay
            if flag
                warning('GMRES returned flag %d',flag)
            end

            % turn x into a wavefunction expansion
            self.soln = unpack(x,orders,pos,self.incident_field.kwave);

            %----------------------------------
            % matrix vector product
            %----------------------------------

            % indented function... this can access all local variables

            function y = matrix_product(x)

                % unpack the vector into wave function expansions
                c = unpack(x,orders,pos,self.incident_field.kwave);

                % loop through destination particles
                for j=1:length(self.particles)

                    % initialise sum to zero
                    csum = regularzero(orders(j),pos(j),self.incident_field.kwave);

                    % loop through source particles
                    for i=1:length(self.particles)

                        % only sum fields from different particles
                        if i~=j

                            % turn the field from particle i into a
                            % regular wavefunction expansion about origin j
                            csum = csum + regularwavefunctionexpansion(c{i},...
                                pos(j),orders(j));

                        end

                    end

                    % temporarily set the origin of the T-matrix so it is
                    % compatible with the wavefunction expansion
                    self.particles(j).tmatrix.setOrigin(pos(j));

                    % apply the T-matrix
                    tsum = self.particles(j).tmatrix * csum;

                    % compute total field scattered by particle j
                    d{j} = c{j} - tsum;

                end

                % convert into a vector
                y = pack(d);

            end

        end

        %----------------------------------
        % farfield function
        %----------------------------------

        function val = getFarfield(self,points)

            % check solution has been computed
            if isempty(self.soln)
                error('First run solve()')
            end

            % compute far field from first particle
            val = self.soln{1}.evaluateFarField(points);

            % add on field from other particles
            for j=2:length(self.particles)
                val = val + self.soln{j}.evaluateFarField(points);
            end

        end

        %----------------------------------
        % scattered field function
        %----------------------------------

        function val = getScatteredField(self,points)

            % check solution has been computed
            if isempty(self.soln)
                error('First run solve()')
            end

            % compute field from first particle
            val = self.soln{1}.evaluate(points);

            % add on field from other particles
            for j=2:length(self.particles)
                val = val + self.soln{j}.evaluate(points);
            end

        end

        %----------------------------------
        % total field function
        %----------------------------------

        function val = getTotalField(self,points)

            % sum scattered and incident fields
            val = self.getScatteredField(points) ...
                + self.incident_field.evaluate(points);

        end

        %----------------------------------
        % add particle
        %----------------------------------

        function addParticle(self,p)

            % check type is correct
            if ~isa(p,'particle')
                error('p must be of class particle');
            end

            % check wavenumber of T-matrix matches incident wave
            if abs(self.incident_field.kwave-p.tmatrix.kwave) > 1e-12
                error('T-matrix wavenumber for p doesn''t match incident field');
            end

            % add particle to list
            self.particles(length(self.particles)+1) = p;
            
        end

        %----------------------------------
        % plot particles
        %----------------------------------

        function varargout = schematic(self)

            % check that some particles exist
            if length(self.particles)==0
                return
            end

            % store hold state so we can restore it later
            hold_state = ishold;

            % plot the first particle
            h = self.particles(1).plot();

            % hold figure for later particles...
            hold on

            % loop through other particles and plot them
            for j=2:length(self.particles)
                tmp = self.particles(j).plot();
                h = [h;tmp];
            end

            % return the graphics handles if required
            if nargout>0
                varargout{1} = h;
            end

        end

        %----------------------------------
        % plot field
        %----------------------------------

        function plot(self,range,phase)

            % set default for phase
            if nargin<3
                phase = 1;
            end               

            % store hold state
            hold_state = ishold;

            % plot the total field... note this may blow up close the the
            % scatterer origins
            tmplot(range,phase,self.incident_field,self.soln{:})

            % next we apply a mask for each scatterer to hide the blow up

            % get the surface plot we just did... it should be the first
            % child of the current axes
            obj = gca().Children(1);

            % check we got what we expected ie a surface
            if ~isa(obj,'matlab.graphics.chart.primitive.Surface')
                error('Couldn''t find expected surface plot')
            end
            
            % extract the x,y and colour data from the surface object
            x = get(obj,'xdata');
            y = get(obj,'ydata');
            c = get(obj,'cdata');

            % write (x,y) in complex format to make next part simpler
            z = x + 1i * y;

            % loop through particles
            for j=1:length(self.particles)
                % if the particle has a radius associated with it
                if isprop(self.particles(j).type,'radius')
                    % mask out parts of c that are inside the radius
                    c(abs(z-self.particles(j).pos)<self.particles(j).type.radius) = NaN;
                end
            end

            % put the masked colour data in the surface plot
            set(obj,'cdata',c)

            % make the plot look pretty
            view([0 90])
            axis equal
            shading interp

            % draw the scatterers too
            hold on
            self.schematic()

            % restore the hold state
            if ~hold_state
                hold off
            end

        end

        %----------------------------------
        % movie of field
        %----------------------------------

        function m = movie(self,range)

            % initialise movie object
            m = mymovie();

            % set number of frames
            n = 26;

            % preliminary run through movie frames to store the axis and 
            % caxis state so we can fix these in the final version
            for f=0:n-1

                % compute the phase
                phase = exp(-1i*f*2*pi/n);

                % plot the figure
                self.plot(range,phase)

                % store axis and caxis
                m.sniff()

            end

            % main run through to generate the frames
            for f=0:n-1

                % compute the phases
                phase = exp(-1i*f*2*pi/n);

                % plot the figure
                self.plot(range,phase)

                % set the axis and xaxis and store the frame
                m.snap()
                
            end

        end

    end

end

%----------------------------------
% pack radiating wavefunction expansion
% coefficients into a vector
%----------------------------------

function vec = pack(c)

% create a temporary cell array
tmp = cell(length(c),1);

% put coefficients into the cell
for j=1:length(c)
    tmp{j} = c{j}.coefficients;
end

% reshape the cell array into a vector
vec = reshape(cell2mat(tmp),[],1);

end


%----------------------------------
% unpack radiating wavefunction expansion
% coefficients from a vector
%----------------------------------

function c = unpack(vec,orders,pos,kwave)

% get index of start of each block in the vector
start = [1;cumsum(2*orders+1)+1];

% put coefficients into radiating wave function expansions
for j=1:length(orders)
    c{j} = radiatingwavefunctionexpansion(orders(j),pos(j),kwave,...
        vec(start(j):start(j+1)-1));
end

end