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

function [d_k_m] = insertion_aux(obj, d_k_m, aux_offsets_f, aux_offsets_t, pilot_pos_f, pilot_pos_t)
    
    % We will apply the mask to a window around the pilot.
    % In order to do so, we need to know the dimensions of this window.
    window_f_size = size(obj.ambiguity_interference,1);
    window_t_size = size(obj.ambiguity_interference,2);
    window_f_offset = (window_f_size-1)/2;
    window_t_offset = (window_t_size-1)/2;
    window_f_center = window_f_offset + 1;
    window_t_center = window_t_offset + 1;
          
    % aux placement has to be done for each pilot
    for i = 1:1:numel(pilot_pos_f)
            
        % extract pilot coordinates within d_k_m
        row = pilot_pos_f(i);
        col = pilot_pos_t(i);        
        
        % from the equalizer equations we need a few values
        k_dot = row - 1;
        m_dot = col - 1;
        k_mat = obj.ambiguity_k_k_dot_diff + k_dot;
        m_mat = obj.ambiguity_m_m_dot_diff + m_dot;
        
        % pre factor
        amb_func_pre_factor = obj.ambiguity_m_m_dot_diff.*k_mat;
        amb_func_pre_factor = (-1).^(mod(amb_func_pre_factor,2));
        
        % pre factor with re im pattern
        amb_func_pre_factor_re_im = amb_func_pre_factor.*(1i).^(mod(k_mat + m_mat, 2));
        
        % determine the exact position of the window within d_k_m
        interference_f_range = row - window_f_offset : row + window_f_offset;
        interference_t_range = col - window_t_offset : col + window_t_offset;
        
        % apply the window for the data
        full_pre_factor = amb_func_pre_factor_re_im.*obj.ambiguity_interference;
        interference_data = full_pre_factor.*d_k_m(interference_f_range,interference_t_range);
        
        % we need to know the sum of all interference
        interference_data = sum(sum(interference_data));
        
        % extract the interference for each auxiliary individually for an assumed amplitude of 1.0
        interference_aux_indiv = zeros(numel(aux_offsets_f), 1);
        for j=1:numel(aux_offsets_f)
            aux_pos_window_f = window_f_center - aux_offsets_f(j);
            aux_pos_window_t = window_t_center - aux_offsets_t(j);
            interference_aux_indiv(j) = full_pre_factor(aux_pos_window_f, aux_pos_window_t)*1.0;
        end
        
        % the actual interference depends on wether our symbol is complex or imaginary
        if mod(k_dot + m_dot, 2) == 0
            % our symbol is real, so the interference will be imaginary
            interference_data = imag(interference_data);
            interference_aux_indiv = imag(interference_aux_indiv);
        else            
            % our symbol is imaginary, so the interference will be real
            interference_data = real(interference_data);
            interference_aux_indiv = real(interference_aux_indiv);
        end
        
        % calculate the weights for each pilot
        weights = equalizer_lib.insertion_aux_solver(interference_data, interference_aux_indiv);
        
        % place each aux with correct weight
        for j=1:numel(aux_offsets_f)
            aux_pos_f = row+aux_offsets_f(j);
            aux_pos_t = col+aux_offsets_t(j);
            d_k_m(aux_pos_f,aux_pos_t) = 1.0*weights(j);
        end
    end
end