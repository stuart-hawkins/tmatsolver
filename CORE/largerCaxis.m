% Determine 'largest' caxis
%
%   cax = largerCaxis(cax1,cax2) computes the smallest caxis cax that contains
%   the caxes cax1 and cax2.
%
%   largerCaxis(cax1,cax2) sets the current caxis to the smallest caxis that 
%   contains the caxes cax1 and cax2.
%
%   cax = largerCaxis(cax1) computes the smallest caxis cax that contains the
%   current caxis and cax1.
%
%   largerCaxis(cax1) sets the current caxis to the smallest caxis that
%   contains the current caxis and cax1.
%
% See also: caxis, largerAxis.
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


function varargout = largerCaxis(varargin)

% if no caxis data is provided then use the current caxis setting
if nargin>1
    cax = varargin{2};
else
    cax = caxis;
end

% put varargin{1} in an ordinary variable so it is easy to access
newcax = varargin{1};

% compare newcax and cax and take the extremal values
newcax(1) = min(newcax(1),cax(1));
newcax(2) = max(newcax(2),cax(2));

% decide what to do with newax...
if nargout==0
    % if no return value then change the current caxis settings
    caxis(newcax)
else
    % return newcax
    varargout{1} = newcax;
end
