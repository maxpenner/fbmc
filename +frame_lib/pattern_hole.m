%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

classdef pattern_hole
    
    properties
        time_start;     % refers to complex frame
        time_end;       % "
        freq_start;     % "
        freq_end;       % "
    end
    
    methods (Static = true, Access = public)
        
        function obj = pattern_hole(time_start_arg, time_end_arg, freq_start_arg, freq_end_arg)
            obj.time_start	= time_start_arg;
            obj.time_end	= time_end_arg;
            obj.freq_start	= freq_start_arg;
            obj.freq_end	= freq_end_arg;
        end
    end
end

