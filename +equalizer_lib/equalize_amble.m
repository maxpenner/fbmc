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

function [amble_ch] = equalize_amble(obj)

    % start with the amble
    if strcmp(obj.amble_pos,'pre')
        
        % determine indices for quick access
        idx_pre = obj.d_k_m_amble_indices.first.lin;
        
        % extract received pilot symbol in preamble
        amb_1_rx = obj.oqam_frame_rx(idx_pre,1);
        amb_2_rx = obj.oqam_frame_rx(idx_pre,3);
        
        % perform zero forcing
        amb_1_zf = amb_1_rx./obj.amble_std(idx_pre);
        amb_2_zf = amb_2_rx./obj.amble_std(idx_pre);
        
        % write to struct
        amble_ch.amb_1_zf = amb_1_zf;
        amble_ch.amb_2_zf = amb_2_zf;

    elseif strcmp(obj.amble_pos,'end')

        % determine indices for quick access
        idx_end = obj.d_k_m_amble_indices.third.lin_relative;
        
        % extract received amble vectors
        amb_3_rx = obj.oqam_frame_rx(idx_end,end-2);
        amb_4_rx = obj.oqam_frame_rx(idx_end,end);
        
        % perform zero forcing
        amb_3_zf = amb_3_rx./(1i*obj.amble_std(idx_end));
        amb_4_zf = amb_4_rx./(1i*obj.amble_std(idx_end));
        
        % write to struct
        amble_ch.amb_3_zf = amb_3_zf;
        amble_ch.amb_4_zf = amb_4_zf;
        
    elseif strcmp(obj.amble_pos,'pre_end')
        
        % determine indices for quick access
        idx_pre = obj.d_k_m_amble_indices.first.lin;
        idx_end = obj.d_k_m_amble_indices.third.lin_relative;
        
        % extract received amble vectors
        amb_1_rx = obj.oqam_frame_rx(idx_pre,1);
        amb_2_rx = obj.oqam_frame_rx(idx_pre,3);
        amb_3_rx = obj.oqam_frame_rx(idx_end,end-2);
        amb_4_rx = obj.oqam_frame_rx(idx_end,end);
        
        % perform zero forcing
        amb_1_zf = amb_1_rx./obj.amble_std(idx_pre);
        amb_2_zf = amb_2_rx./obj.amble_std(idx_pre);
        amb_3_zf = amb_3_rx./(1i*obj.amble_std(idx_end));
        amb_4_zf = amb_4_rx./(1i*obj.amble_std(idx_end));
        
        % write to struct
        amble_ch.amb_1_zf = amb_1_zf;
        amble_ch.amb_2_zf = amb_2_zf;
        amble_ch.amb_3_zf = amb_3_zf;
        amble_ch.amb_4_zf = amb_4_zf;
        
    else
        error('Unknown preamble position.');
    end
end