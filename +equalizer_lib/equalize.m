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

function [oqam_frame_equ_rx, ch_estimation, ch_estimation_inv, warp_fac] = equalize(obj)
    
    % call a specific equalizer
    if strcmp(obj.equalizer_type, 'none')
        
        % all ones means no equalization
        ch_estimation = ones(size(obj.oqam_frame_rx));
        
    elseif strcmp(obj.equalizer_type, 'amble')
        
        % SECURITY
        if obj.amble_len <= 0
            error('Cannot estimate channel with amble if no amble was sent.');
        end
        
        % equalize amble
        amble_ch = equalizer_lib.equalize_amble(obj);
        
        % estimate the entire equalization matrix only from ambles
        ch_estimation = equalizer_lib.equalize_interpolate_amble(obj, amble_ch);

    elseif strcmp(obj.equalizer_type, 'pilot')
                
        % SECURITY
        if numel(obj.pilots) == 0
            error('Cannot estimate channel with pilots if no pilots were sent.');
        end        
        
        % equalize pilots
        pilot_ch = equalizer_lib.equalize_pilot(obj);
        
        % estimate the entire equalization matrix only from pilots
        ch_estimation = equalizer_lib.equalize_interpolate_pilot(obj, pilot_ch);
        
    elseif strcmp(obj.equalizer_type, 'amble_pilot')
        
        % SECURITY
        if obj.amble_len <= 0 || numel(obj.pilots) == 0
            error('Cannot estimate channel with amble and pilots if no amble and/or no pilots were sent.');
        end
        
        % equalize preamble
        amble_ch = equalizer_lib.equalize_amble(obj);
        
        % equalize pilots
        pilot_ch = equalizer_lib.equalize_pilot(obj);
        
        % estimate the entire equalization matrix from ambles and pilots
        ch_estimation = equalizer_lib.equalize_interpolate_amble_pilot(obj, amble_ch, pilot_ch);
        
    elseif strcmp(obj.equalizer_type, 'perfect')
        
        % use the perfect channel knowledge
        ch_estimation = obj.perfect_ch;     
        
    else
        error('Unknown equalizer type.');
    end
    
    % invert the channel
    ch_estimation_inv = 1./ch_estimation;

    % equalize the oqam frame by removing the channel effect
    oqam_frame_equ_rx = obj.oqam_frame_rx.*ch_estimation_inv;
    
    % dewarping has to be done after equalization
    if strcmp(obj.warp_rx,'none') == true
        
        warp_fac = [];
        
    elseif strcmp(obj.warp_rx,'measure') == true
        
        % from the estimated channel we can calculate our warping factors
        warp_fac = warp_lib.warp_rx_measure(obj,ch_estimation);
        
    elseif strcmp(obj.warp_rx,'dewarp') == true
        
        % we have equalized our frame, now we also have to remove the warping that was applied at the transmitter
        oqam_frame_equ_rx = oqam_frame_equ_rx.*(1./obj.warp_fac);
        
        % for now, keep the warp factor
        warp_fac = obj.warp_fac;
    end    
    
    %% Debugging Plot
    if 1==0
        %% plot received oqam frame before equalization
        oqam_frame_rx_abs = pow2db(abs(obj.oqam_frame_rx));

        % limit to a minimum
        minimum = -15;
        maximum = 5;
        oqam_frame_rx_abs(oqam_frame_rx_abs < minimum) = minimum;
        
        figure()
        clf()
        imagesc(oqam_frame_rx_abs)
        axis image
        axis ij
        colormap jet
        colorbar
        xlabel('symbol m dot');
        ylabel('subcarrier index');
        title('Received OQAM frame Absolute before Equalizer');
        
        % adjust scaling only for 2d
        caxis([minimum maximum]);
        
        %% plot channel esimation
        channel_estimation_abs = pow2db(abs(ch_estimation));

        % limit to a minimum
        channel_estimation_abs(channel_estimation_abs < minimum) = minimum;
        
        figure()
        clf()
        imagesc(channel_estimation_abs)
        axis image
        axis ij
        colormap jet
        colorbar
        xlabel('symbol m dot');
        ylabel('subcarrier index');
        title('Measured Channel Absolute from Equalizer');
        
        % adjust scaling only for 2d
        caxis([minimum maximum]);
        
        %% plot equalized frame
        oqam_frame_equ_rx_abs = pow2db(abs(oqam_frame_equ_rx));

        % limit to a minimum
        oqam_frame_equ_rx_abs(oqam_frame_equ_rx_abs < minimum) = minimum;
        
        figure()
        clf()
        imagesc(oqam_frame_equ_rx_abs)
        axis image
        axis ij
        colormap jet
        colorbar
        xlabel('symbol m dot');
        ylabel('subcarrier index');
        title('Equalized OQAM Frame');
        
        % adjust scaling only for 2d
        caxis([minimum maximum]);
        
        %% plot a scatterplot of the equalized symbols
        scatterplot(oqam_frame_equ_rx(obj.d_k_m_data_indices.lin));
        title('Scatterplot of data points');
        
        %% how some numerical results
        fprintf('Average magnitude at data points of channel estimation: %f\n', mean(abs(ch_estimation(obj.d_k_m_data_indices.lin))));
        fprintf('Average magnitude at data points of inverse channel estimation: %f\n', mean(abs(ch_estimation_inv(obj.d_k_m_data_indices.lin))));
    end
end