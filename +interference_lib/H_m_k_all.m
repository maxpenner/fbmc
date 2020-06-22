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

function [H_m_k_cell] = H_m_k_all(obj, channel_handle, mu_max, kappa_max, include_factor)
    
    % we reconstruct the frame before 
    H_m_k_cell = cell(2*kappa_max+1, 2*mu_max+1);
    
    % now sum up the interference
    idx_f = 1;
    for kappa = -kappa_max : 1 : kappa_max
        
        idx_t = 1;
        
        % each time offset relevant
        for mu = -mu_max : 1 : mu_max
            
            H_m_k_cell(idx_f,idx_t) = {interference_lib.H_m_k(obj, channel_handle, mu, kappa, include_factor)};
            
            idx_t = idx_t + 1;
        end
        
        idx_f = idx_f + 1;
    end
    
%     % check if frames are equal
%     if strcmp(include_factor,'include') == true
%         % do nothing
%     elseif strcmp(include_factor,'include not') == true
%         error('Cannot test oqam frame reconstruction if factor in H m k is ommited.');
%     else
%         error('Unknown factor behaviour.');
%     end
%     
%     % now reconstruct the actual frame
%     oqam_frame_reconstruct_before_equ = zeros(size(obj.oqam_frame));
% 
% 	% now sum up the interference
%     idx_1 = 1;
%     for kappa = -kappa_max : 1 : kappa_max
%         
%         idx_0 = 1;
%         
%         % each time offset relevant
%         for mu = -mu_max : 1 : mu_max
%             
%             % time shift
%             if mu >= 0
%                 oqam_frame_shifted = [obj.oqam_frame(:,1+mu:end), zeros(obj.n_subc,mu)];
%             else
%                 oqam_frame_shifted = [zeros(obj.n_subc,-mu), obj.oqam_frame(:,1:end+mu)];
%             end
%             
%             % freq shift
%             oqam_frame_shifted = circshift(oqam_frame_shifted,-kappa,1);
%             
%             % extract correct factors
%             H_m_k = H_m_k_cell{idx_1,idx_0};
%             
%             % reconstruct the oqam frame
%             oqam_frame_reconstruct_before_equ = oqam_frame_reconstruct_before_equ + oqam_frame_shifted.*H_m_k;
%             
%             idx_0 = idx_0 + 1;
%         end
%         idx_1 = idx_1 + 1;
%     end
%     
%     error_abs = abs(obj.oqam_frame_rx - oqam_frame_reconstruct_before_equ);
%     error_max = max(max(error_abs));
%     if error_max > 1e-3
%         error('Large offset.');
%     end
end

