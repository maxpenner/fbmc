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

function [signal_power, error_power, sinr_dB, bins_counts, bin_centers] = points_real_imag_sinr_hist_db(obj, point_type, n_bins, edge_l, edge_r)

    % SECURITY
    if strcmp(obj.modulation_type, 'PAM') == false
        error('SINR on real points can only be measured is we are using PAM');
    end        

    % extract correct real points from both d_k_m versions
    if strcmp(point_type, 'data')
        
        points_tx = obj.d_k_m(obj.d_k_m_data_indices.lin);
        points_rx = obj.d_k_m_rx(obj.d_k_m_data_indices.lin);
        
        histogramm_vector = points_rx - points_tx;
        
    elseif strcmp(point_type, 'pilot')
        
        points_tx = obj.d_k_m(obj.d_k_m_pilot_indices.lin);
        points_rx = obj.d_k_m_rx(obj.d_k_m_pilot_indices.lin);
        
        histogramm_vector = points_rx - points_tx;
        
    elseif strcmp(point_type,'data_oqam_frame_equ_rx')
        
        points_tx = obj.oqam_frame(obj.d_k_m_data_indices.lin);
        points_rx = obj.oqam_frame_equ_rx(obj.d_k_m_data_indices.lin);
        
        histogramm_vector = abs(points_rx - points_tx);
        
    elseif strcmp(point_type, 'data_warped')
        
        % At the transmitter, the data points where warped, see frame_comp.m.
        % But we have to compare to an unwarped version of those data points.
        % This is possible because we assume both the transmitter and the receiver know the warping factor.
        d_k_m_unwarped = obj.d_k_m./obj.warp_fac;
        
        points_tx = d_k_m_unwarped(obj.d_k_m_data_indices.lin);
        points_rx = obj.d_k_m_rx(obj.d_k_m_data_indices.lin);
        
        histogramm_vector = points_rx - points_tx;
        
    else
        error('Unknown real point type.');
    end
    
    % get power of transmitted points = rms^2 = Effektivwert^2
	signal_power = rms(reshape(points_tx,1,[]))^2;
    
    % determine euclidian distances between send real points and received real points
    error_euclid_distance = abs(points_rx - points_tx);
        
    % get error power = rms^2 = Effektivwert^2
	error_power = rms(reshape(error_euclid_distance,1,[]))^2;

	% convert to dB
   	sinr_dB = pow2db(signal_power/error_power);
    
    % SECURITY
    if edge_r <= edge_l
        error('Left edge must be smaller than right edge.');
    end

    % first calculate distance between edges
    [edges, bin_centers, ~] = interference_lib.hist_edges(n_bins, edge_l, edge_r);
    
    % reshape to vector
    histogramm_vector = reshape(histogramm_vector, 1, []);
    
    % let matlab perform the histogramm
    [bins_counts, ~] = histcounts(histogramm_vector, edges);
end

