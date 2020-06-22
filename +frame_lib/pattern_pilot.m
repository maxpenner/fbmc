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

%
%   QAM:
%
%       DD      AX
%       PA      PX
%       DD      AX
%
%
%   PAM:
%
%     DDAD
%     DAPA      any of the four can be active of inactive
%     DDAD
%
%   Pilot indices:
%
%       1
%      4P2
%       3

classdef pattern_pilot
    
    properties
        
        % position of the actual pilot
        time_offset;    % refers to complex frame
        freq_offset     % "
        time_step;      % "
        freq_step;      % "
        
        % which constellation are we using
        constellation_type;
        
        % index of each pilot
        aux_indices;
        
        % offsets relative to pilot of the actual auxiliary pilots
        t_aux_offset;
        f_aux_offset;
        
        % We can block certain complex points that are neither data nor pilot.
        % Only relevant is we use QAM instead of PAM.
        block_t_offset;
        block_f_offset;
    end
    
    methods (Static = true, Access = public)
        
        function obj = pattern_pilot(time_offset_arg, freq_offset_arg, time_step_arg, freq_step_arg, constellation_type_arg, aux_indices_arg)
            
            obj.time_offset         = time_offset_arg;
            obj.freq_offset         = freq_offset_arg;
            obj.time_step           = time_step_arg;
            obj.freq_step           = freq_step_arg;
            obj.constellation_type  = constellation_type_arg;            
            obj.aux_indices         = aux_indices_arg;
            
            % make sure pilots are not doubling
            if numel(obj.aux_indices) ~= numel(unique(obj.aux_indices))
                error('Pilots cannot be inserted twice.');
            end
            
            % make sure pilots are ascending
            if obj.aux_indices ~= sort(obj.aux_indices)
                error('Pilots must be given in ascending order.');
            end            
            
            % if we are using QAM, we always have to use a pair of real and imaginary part
            if strcmp(obj.constellation_type, 'QAM')
                
                % check if combination of pilots is possible
                if numel(obj.aux_indices) == 0
                    
                    error('We must use pilots 2 or 1 and 3.');
                    
                elseif numel(obj.aux_indices) == 1
                    
                    if obj.aux_indices ~= 2
                        error('With one pilot only pilot 2 is possible.');
                    else
                        obj.t_aux_offset = 1;
                        obj.f_aux_offset = 0;
                        obj.block_t_offset = [];
                        obj.block_f_offset = [];
                    end
                    
                elseif numel(obj.aux_indices) == 2
                    
                    if obj.aux_indices(1) ~= 1 || obj.aux_indices(2) ~= 3
                        error('With two pilots pilots only the combination 1 and 3 is possible. Also pilots have to be sorted.');
                    else
                        obj.t_aux_offset = [0, 0];
                        obj.f_aux_offset = [-1, 1];
                        obj.block_t_offset = [1, 1, 1];
                        obj.block_f_offset = [-1, 0, 1];                        
                    end
                    
                else
                    error('Incorrect number of pilots.');
                end
                
            % if we are using PAM, any re and im are independent
            elseif strcmp(obj.constellation_type, 'PAM')
                
                % check if combination of pilots is possible
                if numel(obj.aux_indices) == 0
                    
                    error('We must use at least aux pilot 1, 2, 3 or 4.');
                    
                elseif numel(obj.aux_indices) >= 1 && numel(obj.aux_indices) <= 4
                    
                    % transform aux indices to an array of offsets
                    obj.t_aux_offset = zeros(numel(obj.aux_indices),1);
                    obj.f_aux_offset = zeros(numel(obj.aux_indices),1);
                    idx = 0;
                    for aux_idx = obj.aux_indices
                        idx = idx + 1;
                        switch aux_idx
                            case 1
                                obj.t_aux_offset(idx) = 0;
                                obj.f_aux_offset(idx) = -1;
                            case 2
                                obj.t_aux_offset(idx) = 1;
                                obj.f_aux_offset(idx) = 0;
                            case 3
                                obj.t_aux_offset(idx) = 0;
                                obj.f_aux_offset(idx) = 1;
                            case 4
                                obj.t_aux_offset(idx) = -1;
                                obj.f_aux_offset(idx) = 0;                                
                            otherwise
                                error('Only aux positions 1, 2, 3 and 4 are allowed.');
                        end
                    end 
                    
                    obj.block_t_offset = [];
                    obj.block_f_offset = [];
                    
                else
                    error('Incorrect number of pilots.');
                end                
            end
        end
    end
end
