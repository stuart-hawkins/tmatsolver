% Plot wave field
%
%  tmplot(ax,[],u) plots the field u over the rectangular domain 
%  [ax(1) ax(2)] x [ax(3) ax(4)]. Here u must be of class incident or
%  wavefunctionexpansion. The quantity plotted is real(u).
%
%  tmplot(ax,[],u,v,...) plots the sum u + v + ....
%
%  tmplot(ax,[],phase,u,v,...) where phase is a complex scalar plots 
%  phase * (u + v + ...).
%
%  tmplot(ax,mode,...) specifies the kind of plot. mode can be 'ABS' or
%  'REAL' or 'LOG' or a function handle. If mode is a function handle then
%  the result of the function is plotted.
%
%  h = tmplot(...) returns the graphics handle h of the plot.
%
% See also: plane_wave, regularwavefunctionexpansion,
% radiatingwavefunctionexpansion.
%
% Stuart C. Hawkins - 18 March 2023

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


function varargout = tmplot(range,mode,varargin)

%-----------------------------------------------
% process the parameters
%-----------------------------------------------

% check range has the right size
if length(range) ~= 4
    error('range must be a vector of length 4')
end

% set default for mode
if isempty(mode)
    mode = 'REAL';
end

% check the phase
if length(varargin)>0 && isnumeric(varargin{1})
    % then varargin{1} is the phase
    phase = varargin{1};

    % wave-arguments are in varargin{2:end}
    varargin_offset = 2;
else
    % default phase is 1
    phase = 1;

    % wave-arguments are in varargin{1:end}
    varargin_offset = 1;
end

% check we have some wave-type parameters passed so we can work out the
% wavenumber
if length(varargin)<varargin_offset
    error('At least one argument of class incident or wavefunctionexpansion must be provided')
end

%-----------------------------------------------
% set constants
%-----------------------------------------------

% points per wavelength parameter
ppw = 20;

%-----------------------------------------------
% try to determine wavenumber
%-----------------------------------------------

% We need the wavenumber because we base the resolution
% of the plot so that we have ppw points per wavelength

% initialise wavenumber
kwave = 0;

% try to determine wavenumber from the fields we have been given
% Note: this bit of code will break down if we have only an incident field
% that is the sum of many objects of class incident
for j=varargin_offset:length(varargin)
    if isa(varargin{j},'wavefunctionexpansion')
        kwave = varargin{j}.kwave;
        break
    elseif isa(varargin{j},'plane_wave')
        kwave = varargin{j}.kwave;
        break
    elseif isa(varargin{j},'point_source')
        kwave = varargin{j}.kwave;
        break
    elseif isa(varargin{j},'incidentplus') || ...
            isa(varargin{j},'incidentminus') || ...
            isa(varargin{j},'incidenttimes')
        kwave = varargin{j}.left.kwave;
        break
    end
end

% if kwave is still zero then we know we weren't able to get it from the
% parameters passed... raise an error
if kwave==0
    error('Unable to determine wavenumber from provided arguments')
end

% get wavelength
lambda = 2*pi/kwave;

%-----------------------------------------------
%  setup grid for plotting
%-----------------------------------------------

% get size of domain in wavelengths
sx = (range(2)-range(1))/lambda;
sy = (range(4)-range(3))/lambda;

% work out how many points to use in each direction
nx = floor(ppw * sx);
ny = floor(ppw * sy);

% ensure minimum number of points (necessary if the wavelength is very
% long)
nx = max(nx,40);
ny = max(ny,40);

% check number of points isn't too large
if nx > 400
    warning('Requested %d points in the x-direction but capped at 400',nx)
    nx = 400;
end
if ny > 400
    warning('Requested %d points in the y-direction but capped at 400',ny)
    ny = 400;
end

% setup grid on which to plot the field
u = linspace(range(1),range(2),nx);
v = linspace(range(3),range(4),ny);
[x,y] = meshgrid(u,v);
z = x + 1i * y;

%-----------------------------------------------
% compute the field and plot
%-----------------------------------------------

% compute field
% Note: we take advantage that wavefunctionexpansion and incident classes
% both have an evaluate method
u = phase * varargin{varargin_offset}.evaluate(z);
for j=varargin_offset+1:length(varargin)
    u = u + phase * varargin{j}.evaluate(z);
end

% compute the quantity to be plotted... what to do depends on mode
if isa(mode,'function_handle')
    % mode is a function that we apply
    v = mode(u);
elseif isa(mode,'char')
    % mode is a string specifying what to do
    switch upper(mode)
        case 'ABS'
            v = abs(u);
        case 'REAL'
            v = real(u);
        case 'IMAG'
            v = imag(u);
        case 'LOG'
            v = log10(abs(u));
        otherwise
            error('mode %s not recognised',mode)
    end
else
    % mode is not the correct type
    error('mode must be a string or a function handle')
end

%-----------------------------------------------
% plot the field
%-----------------------------------------------

% plot field
h = surf(x,y,zeros(size(u)),v);

% return the graphics handle if required
if nargout>0
    varargout{1} = h;
end
