% Determine 'largest' axis
%
%   ax = largerAxis(ax1,ax2) computes the smallest axis ax that contains
%   the axes ax1 and ax2.
%
%   largerAxis(ax1,ax2) sets the current axis to the smallest axis that 
%   contains the axes ax1 and ax2.
%
%   ax = largerAxis(ax1) computes the smallest axis ax that contains the
%   current axis and ax1.
%
%   largerAxis(ax1) sets the current axis to the smallest axis that
%   contains the current axis and ax1.
%
% See also: axis, largerCaxis.
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


function varargout = largerAxis(varargin)

% if no axis data is provided then use the current axis setting
if nargin>2
    ax = varargin{2};
else
    ax = axis;
end

% put varargin{1} in an ordinary variable so it is easy to access
newax = varargin{1};

% compare newax and ax and take the extremal values
newax(1) = min(newax(1),ax(1));
newax(3) = min(newax(3),ax(3));
newax(2) = max(newax(2),ax(2));
newax(4) = max(newax(4),ax(4));

% decide what to do with newax...
if nargout==0
    % if no return value then change the current axis settings
    axis(newax)
else
    % return newax
    varargout{1} = newax;
end
