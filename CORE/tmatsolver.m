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
%  a = obj.getCoefficients() returns the vector of coefficients for the
%  wavefunction expansions.
%
%  a = obj.getCoefficients(k) returns the vector of coefficients for the
%  kth scatterer.
%
%  a = obj.getCoefficients(l,k) returns the order l coefficient for the
%  kth scatterer.
%
%  obj.schematic() visualises the scatterers.
%
%  obj.plot(ax) plots the real part of the total field on the domain
%  [ax(1) ax(2)] x [ax(3) ax(4)].
%
%  obj.plot(ax,phase) plots real(phase*u) where u is the total field.
%
%  obj.plot(ax,phase,'ABS') plots abs(phase*u) where u is the total field.
%
%  obj.plot(ax,phase,'LOG') plots log10(abs(phase*u)) where u is the total 
%  field.
%
%  obj.plot(ax,phase,f) plots f(phase*u) where u is the total field and f
%  is a function handle.
%
%  obj.plotIncident(...) plots the incident field.
%
%  obj.plotScattered(...) plots the scattered field.
%
%  m = obj.movie(ax) generate a movie of the real part of the total field 
%  on the domain [ax(1) ax(2)] x [ax(3) ax(4)]. The movie is stored in an 
%  object of type mymovie.
%
% The solve() method uses GMRES to solve a linear system. Additional
% methods are provided to control and monitor the performance of the
% solver.
%
%  obj.setSolverTol(tol) sets the solver tolerance to tol.
%
%  obj.setSolverIterations(N) sets the maximum number of iterations to N.
%
%  obj.getSolverIterations() returns the number of iterations that GMRES
%  required to match the solver tolerance.
%
%  obj.getSolverResidualHistory() returns a vector containing the GMRES
%  residual norm at each iteration.
%
%  obj.getSolverFlag() returns the GMRES flag. See help gmres for help on
%  interpreting this.
%
% See also: particle, plane_wave, mymovie.
%
% Stuart C. Hawkins - 18 March 2023

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
        solver_tol
        solver_itns
        solver_actual_iterations
        solver_reshist
        solver_flag
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

            % set default GMRES tolerance
            self.solver_tol = 1e-8;

            % set default GMRES restart parameter
            self.solver_itns = [];

            % initalise residual history, flag and iteration number
            self.solver_actual_iterations = [];
            self.solver_flag = [];
            self.solver_reshist = [];

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
            dof = sum(2*orders+1);

            % initialise vector holding centres of the particles
            pos = zeros(length(self.particles),1);

            % initialise vector holding rotation of the particles
            rot = zeros(length(self.particles),1);

            % pull out centres and rotations
            for j=1:length(self.particles)
                pos(j) = self.particles(j).pos;
                rot(j) = self.particles(j).rot;
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

                % rotate the coordinate system to match the scatterer
                a{j}.rotatecoordinates(rot(j))

                % apply the T-matrix to the incident field
                b{j} = self.particles(j).tmatrix * a{j};

                % rotate the coordinate system back again
                b{j}.rotatecoordinates(-rot(j))

            end

            % pack the coefficients of b into a vector
            rhs = pack(b);

            %- - - - - - - - - - - - - - - - -
            % setup solver iterations parameters
            %- - - - - - - - - - - - - - - - -

            % we choose the GMRES parameters to try to avoid it restarting

            if isempty(self.solver_itns)

                % set nitns to empty... this means the number of iterations
                % is set by maxit
                nitns = [];

                % set maxit to the maximum possible ie dof
                maxit = dof;

            else

                % set number of GMRES iterations
                nitns = self.solver_itns;

                % and no restart
                maxit = 1;

            end

            %- - - - - - - - - - - - - - - - -
            % solve linear system
            %- - - - - - - - - - - - - - - - -

            % solve linear system using GMRES
            [x,self.solver_flag,relres,self.solver_actual_iterations,self.solver_reshist] = gmres(@matrix_product,rhs,...
                nitns,self.solver_tol,maxit);

            % check GMRES worked okay
            if self.solver_flag
                warning('GMRES returned flag %d',self.solver_flag)
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

                    % rotate the coordinate system to match the scatterer
                    csum.rotatecoordinates(rot(j))
                
                    % apply the T-matrix
                    tsum = self.particles(j).tmatrix * csum;

                    % rotate the coordinate system back again
                    tsum.rotatecoordinates(-rot(j))

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

        function plot(self,range,phase,mode)

            % set default for phase
            if nargin<3 || isempty(phase)
                phase = 1;
            end               

            % set default for mode
            if nargin<4 || isempty(mode)
                mode = 'REAL';
            end

            % store hold state
            hold_state = ishold;

            % plot the total field... note this may blow up close the the
            % scatterer origins
            tmplot(range,mode,phase,self.incident_field,self.soln{:})

            % next we apply a mask for each scatterer to hide the blow up
            self.plotMasks()

            % make the plot look pretty
            view([0 90])
            axis equal
            shading interp

            % draw the scatterers too
            hold on
            self.schematic()

            % make axis tight
            axis(range)

            % restore the hold state
            if ~hold_state
                hold off
            end

        end

        %----------------------------------
        % plot scattered field
        %----------------------------------

        function plotScattered(self,range,phase,mode)

            % set default for phase
            if nargin<3 || isempty(phase)
                phase = 1;
            end               

            % set default for mode
            if nargin<4 || isempty(mode)
                mode = 'REAL';
            end

            % store hold state
            hold_state = ishold;

            % plot the total field... note this may blow up close the the
            % scatterer origins
            tmplot(range,mode,phase,self.soln{:})

            % next we apply a mask for each scatterer to hide the blow up
            self.plotMasks()

            % make the plot look pretty
            view([0 90])
            axis equal
            shading interp

            % draw the scatterers too
            hold on
            self.schematic()

            % make axis tight
            axis(range)

            % restore the hold state
            if ~hold_state
                hold off
            end

        end

        %----------------------------------
        % plot scattered field
        %----------------------------------

        function plotIncident(self,range,phase,mode)

            % set default for phase
            if nargin<3 || isempty(phase)
                phase = 1;
            end               

            % set default for mode
            if nargin<4 || isempty(mode)
                mode = 'REAL';
            end

            % store hold state
            hold_state = ishold;

            % plot the total field... note this may blow up close the the
            % scatterer origins
            tmplot(range,mode,phase,self.incident_field)

            % next we apply a mask for each scatterer to hide the blow up
            self.plotMasks()

            % make the plot look pretty
            view([0 90])
            axis equal
            shading interp

            % draw the scatterers too
            hold on
            self.schematic()

            % make axis tight
            axis(range)

            % restore the hold state
            if ~hold_state
                hold off
            end

        end

        %----------------------------------
        % plot scatterer masks
        %----------------------------------

        % Note: this assumes we just did a plot... other wise strange 
        % results will occur

        function plotMasks(self)

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

        %----------------------------------
        % get coefficients
        %----------------------------------

        function varargout = getCoefficients(self,varargin)

            % check solution has been computed
            if isempty(self.soln)
                error('First run solve()')
            end

            switch nargin

                case 1
                    % then return all the coefficients as a cell array
                    varargout{1} = pack(self.soln);

                case 2
                    % then varargin{1} represents a scatterer's index and
                    % return the corresponding coefficients
                    varargout{1} = pack(self.soln(varargin{1}));

                case 3
                    % then varargin{2} represents a scatterer's index and
                    % varargin{1} the order of the desired coefficient
                    tmp = self.soln{varargin{2}};
                    varargout{1} = tmp.coefficients(varargin{1} + 1 + ...
                        self.particles(varargin{2}).tmatrix.order);

                otherwise

                    error('Unexpected number of parameters')

            
            end

        end

        %----------------------------------
        % get GMRES residual
        %----------------------------------

        function val = getSolverResidualHistory(self)

            % check that the residual history exists
            if isempty(self.solver_reshist)
                error('First call solve method.')
            end

            % return residual history
            val = self.solver_reshist;

        end

        %----------------------------------
        % get GMRES iterations
        %----------------------------------

        function val = getSolverIterations(self)

            % check that the itns exists
            if isempty(self.solver_actual_iterations)
                error('First call solve method.')
            end

            % return solver iterations
            val = self.solver_actual_iterations(end);

        end

        %----------------------------------
        % get GMRES flag
        %----------------------------------

        function val = getSolverFlag(self)

            % check that the flag exists
            if isempty(self.solver_flag)
                error('First call solve method.')
            end

            % return solver flag
            val = self.solver_flag;

        end

        %----------------------------------
        % set GMRES iterations
        %----------------------------------

        function setSolverIterations(self,itns)

            % set default value
            if nargin<2
                itns = [];
            end

            % set solver iterations
            self.solver_itns = itns;

        end

        %----------------------------------
        % set GMRES tolerance
        %----------------------------------

        function setSolverTol(self,tol)

            % set default value
            if nargin<2
                tol = 1e-8;
            end

            % set solver tolerance
            self.solver_tol = tol;

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