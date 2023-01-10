% Class for recording and exporting movies
%
%   The main purpose of this class is to avoid the axis scale and colour
%   axis scale (caxis) changing during a movie. This is avoided by
%   generating the movie twice... the first time so that an axis scale and
%   caxis scale can be fixed, the second time to recorde the frames with
%   the new axis and caxis data.
%
%   m = mymovie() initialises the object.
%
%   m.sniff() takes a snapshot of the axis and caxis properties.
%
%   m.snap() applies the fixed axis and caxis scale and snaps a movie
%   frame.
%
%   m.play() plays back the movie
%
%   m.export(fname,'mp4') exports the movie in mp4 format to filename.
%
%   m.export(fname,'avi') exports the movie in AVI format to filename.
%
%   m.export(fname,'gif') exports the movie gif format to filename.
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


classdef mymovie < handle
    
    properties
        frames
        fps
        length
        ax
        cax
    end
    
    methods
   
        %-------------------------------------------------
        % constructor
        %-------------------------------------------------

        function self = mymovie()
            
            self.frames = struct('cdata',[],'colormap',[]);
            self.length = 0;
            self.fps = 20;
                        
        end
        
        %-------------------------------------------------
        % detect snapshot properties like axis and caxis
        %-------------------------------------------------

        function sniff(self)
                        
            % record axis data
            if isempty(self.ax)
                % first time... store axis
                self.ax = axis;
            else
                % store larger of current stored version and current axis
                self.ax = largerAxis(self.ax);
            end            
            
            % record caxis data
            if isempty(self.cax)
                % first time... store caxis
                self.cax = caxis;
            else
                % store large of current stored version and current caxis
                self.cax = largerCaxis(self.cax);
            end
            
            drawnow

        end
        
        %-------------------------------------------------
        % apply stored axis and caxis and then snap the frame
        %-------------------------------------------------

        function snap(self)
                        
            % apply the stored axis
            if ~isempty(self.ax)
                axis(self.ax)
            end

            % apply the stored caxis
            if ~isempty(self.cax)
                caxis(self.cax)
            end

            % drawnow to make sure the changed axis stuff is applied
            drawnow

            % store the frame
            self.frames(self.length+1) = getframe(gcf);
            self.length = self.length+1;
            
        end
       
        %-------------------------------------------------
        % set movie frames per second
        %-------------------------------------------------

        function set.fps(self,fps)
            
            self.fps = fps;
            
        end
        
        %-------------------------------------------------
        % playback the stored movie
        %-------------------------------------------------

        function play(self,repeat)
           
            % set default repeat value
            if nargin < 2
                repeat = 1;
            end
            
            % make the current figure visible
            figure(gcf)
            
            % plot the movie
            movie(self.frames,repeat,self.fps)
            
        end
        
        %-------------------------------------------------
        % export the stored movie
        %-------------------------------------------------

        function export(self,filename,opts)
            
            % set default format
            if nargin<3 || isempty(opts)
                opts = 'avi';
            end
            
            % output in AVI format
            if contains(opts,'avi')            
                movie2avi(self.frames,filename,'fps',self.fps)                
            end
                
            % output in mp4 format
            if contains(opts,'mp4')
                tmp = VideoWriter(filename,'MPEG-4');
                set(tmp,'Quality',100);
                tmp.FrameRate = self.fps;
                tmp.open;
                tmp.writeVideo(self.frames);
                tmp.close;
            end

            % output as a GIF animation
            if contains(opts,'gif')
                
                for j = 1:self.length
                    [A,map] = rgb2ind(frame2im(self.frames(j)),256);
                    if j == 1
                        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1/self.fps);
                    else
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1/self.fps);
                    end
                end
            end
            
        end
        
    end
    
end