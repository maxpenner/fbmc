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

function [perfect_ch] = perfect_ch_knowledge(obj, channel_handle)

%     % needed when testing sim_model_analytical_A_formula_mod_demod.m
%     perfect_ch = [];
%     return;

    % determine ideal channel knowledge
    if strcmp(channel_handle.type, 'awgn')

        perfect_ch = ones(size(obj.oqam_frame));

    elseif strcmp(channel_handle.type, 'deterministic')

        perfect_ch = ones(size(obj.oqam_frame));
        
        warning('ON');
        warning('Calculating perfect channel no implemented for deterministic channel.');
        warning('ON');

    elseif strcmp(channel_handle.type, 'statistical')

        % first make sure we actually have path gains
        if numel(channel_handle.s_rayleigh_gains) == 0
            error('Cannot calculate H_k_dot_m_dot for statistical channel without path gains.');
        end
        
        % estimate perfect channel
        perfect_ch = interference_lib.H_m_k(obj, channel_handle, 0, 0, 'include');

    else
        error('Unknown channel type.');
    end
end