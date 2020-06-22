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
%   We are at position k_dot and m_dot.
%   The matrix 'interference_matrix' shall contain the interference that neighbors cause at our position.
%   The neighbors are at position k and m.
%   The matrices 'k_k_dot_diffs' and 'm_m_dot_diffs' contain the neighbor position relative to us in the center.
%
%       Arguments:
%           k_k_dot_diff_max	number of subcarrier frequencies
%           m_m_dot_diff_max	number of offsets in multiples of n_subc/2
%           f_res             	number of points within one subcarrier spacing
%           t_res            	number of points within a range of n_subc/2
function [interference_matrix, ...
            k_k_dot_diffs, ...
            m_m_dot_diffs] = ambiguity_interference(obj, k_k_dot_diff_max, m_m_dot_diff_max, f_res, t_res)

    % we create a scaled version of the filter function so that we reach a higher integration precision
    upscale = 1;
    n_subc_scaled = obj.n_subc*upscale;
    
%     % limit the upscaling, we don't need endless precision
%     if n_subc_scaled >= 1024
%         n_subc_scaled = 1024;
%     end
    
    % Re-extract the exact same filter but with better resolution.
    % In the equation the filter_func is p[t].
    filter_func = fbmc_lib.filter_func_load(obj.prototype_filter, obj.n_subc, obj.k);

    % determine all frequencies offsets of interest
	k_k_dot_diffs = -k_k_dot_diff_max : 1 : k_k_dot_diff_max;
    
    % overwrite frequency offsets for higher resolution
    if f_res > 1
        
        % we have the same k_k_dot range, but with higher resolution
        f_offsets_l = linspace(-k_k_dot_diff_max, 0, f_res*k_k_dot_diff_max);
        f_offsets_r = linspace(0, k_k_dot_diff_max, f_res*k_k_dot_diff_max);
        
        % drop the second 0 and combine
        k_k_dot_diffs = [f_offsets_l f_offsets_r(2:end)];
        
        % make sure we still have the zero in the center
        if k_k_dot_diffs((numel(k_k_dot_diffs)-1)/2 + 1) ~= 0
            error('For higher resolution there is no zero in the center of the frequencies array.');
        end
    end
    
    % transpose the frequency, makes more sense if they go from top to bottom
    k_k_dot_diffs = k_k_dot_diffs';
    
    % in the equations we need K, the number of subcarriers
    K_scaled = n_subc_scaled;
    
    % calculate actual mixing frequencies from the k_k_dot differences
    actual_freqs = k_k_dot_diffs*1/K_scaled;

    % create an array for convolution
    filter_func_mixed = repmat(filter_func, numel(actual_freqs), 1);
    
    % Create a time base starting at 0.
    % In the formulas this value is x.
    x = 0:1:numel(filter_func);
    x = x(1:end-1);
    
    % mix each filter function with it's corresponding frequency
    for i=1:1:numel(actual_freqs)
        filter_func_mixed(i,:) = filter_func.*exp(1i*2*pi*actual_freqs(i)*x);
    end

    % now convolute with initial filter function, by that we calculate each possible time offsets
    filter_func_mixed_conv = conv2(filter_func_mixed, fliplr(filter_func));

%     % normalize to maximum in the center
%     normalization_fac = max(max(filter_func_mixed_conv));
%     filter_func_mixed_conv = filter_func_mixed_conv/normalization_fac;

    % determine the center of the array in the horizontal time
    conv_center_index = (size(filter_func_mixed_conv,2)-1)/2+1;

    % determine all points in time of interest
    m_m_dot_diffs = -m_m_dot_diff_max : 1 : m_m_dot_diff_max;    
    
    % overwrite if we want a higher resolution
    if t_res > 1
        
        % we have to make sure that t_resoltion is a divider of the number of subcarriers
        if mod(n_subc_scaled/2, t_res) ~= 0
            error('Resolution of time is not possible. It has to be a divider of scaled n_subc.');
        end
        
        % Overwrite the offsets of interest.
        % Not really a nice solution due to the none integer division, but I hope it works.
        m_m_dot_diffs = -m_m_dot_diff_max : 1/t_res : m_m_dot_diff_max;
    end
    
    % in the formulas we need a factor T
    T_scaled = n_subc_scaled;    
    
    % actual time differences we will extract
    actual_time_offsets = m_m_dot_diffs*T_scaled/2;    
    
    % we have to consider the indices from the center
    t_array_indices = actual_time_offsets + conv_center_index;

    % extract those frequency values we are interested in
    interference_matrix = filter_func_mixed_conv(:, t_array_indices);
    
    % create matrices of k_k_dot and m_m_dot that have the same size as interference_matrix
    k_k_dot_diffs = repmat(k_k_dot_diffs, 1, size(interference_matrix,2));
    m_m_dot_diffs = repmat(m_m_dot_diffs, size(interference_matrix,1), 1);
end
