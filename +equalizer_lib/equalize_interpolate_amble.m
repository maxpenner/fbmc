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

function [channel_estimation] = equalize_interpolate_amble(obj, amble_ch)

    % channel estimation depends on our frame structure
    if strcmp(obj.amble_pos,'pre')
        
        % determine indices for quick access
        idx_pre = obj.d_k_m_amble_indices.first.lin;        
        
        % create an average estimation of both vector
        avg_first_second = (amble_ch.amb_1_zf + amble_ch.amb_2_zf)/2;
        
        % interpolate real and imag part
        if strcmp(obj.interpol_ri_ma, 'ri')
            
            % interpolate over all subcarriers
            amble_channel_estimation_interpolated = interp1(idx_pre, avg_first_second, (1:1:obj.n_subc)', 'spline');

            % create the channel equalization matrix
            channel_estimation = repmat(amble_channel_estimation_interpolated(:,1),1,obj.n_symbols*2);            
            
        elseif strcmp(obj.interpol_ri_ma, 'ma')
            
            % interpolate over all subcarriers
            amble_channel_estimation_interpolated_magnitude = interp1(idx_pre, abs(avg_first_second), (1:1:obj.n_subc)', 'spline');
            amble_channel_estimation_interpolated_angle = interp1(idx_pre, avg_first_second, (1:1:obj.n_subc)', 'spline');            
            amble_channel_estimation_interpolated_angle = angle(amble_channel_estimation_interpolated_angle);

            % combine both
            amble_channel_estimation_interpolated = amble_channel_estimation_interpolated_magnitude.*exp(1i*amble_channel_estimation_interpolated_angle);            

            % create the channel equalization matrix
            channel_estimation = repmat(amble_channel_estimation_interpolated(:,1),1,obj.n_symbols*2);          
            
        end

    elseif strcmp(obj.amble_pos,'end')

        % determine indices for quick access
        idx_end = obj.d_k_m_amble_indices.third.lin_relative;
        
        % create an average estimation of both vector
        avg_third_fourth = (amble_ch.amb_3_zf + amble_ch.amb_4_zf)/2;
        
        % interpolate real and imag part
        if strcmp(obj.interpol_ri_ma, 'ri')
            
            % interpolate over all subcarriers
            amble_channel_estimation_interpolated = interp1(idx_end, avg_third_fourth, (1:1:obj.n_subc)', 'spline');

            % create the channel equalization matrix
            channel_estimation = repmat(amble_channel_estimation_interpolated(:,1),1,obj.n_symbols*2);            
            
        elseif strcmp(obj.interpol_ri_ma, 'ma')
            
            % interpolate over all subcarriers
            amble_channel_estimation_interpolated_magnitude = interp1(idx_end, abs(avg_third_fourth), (1:1:obj.n_subc)', 'spline');
            amble_channel_estimation_interpolated_angle = interp1(idx_end, avg_third_fourth, (1:1:obj.n_subc)', 'spline');            
            amble_channel_estimation_interpolated_angle = angle(amble_channel_estimation_interpolated_angle);

            % combine both
            amble_channel_estimation_interpolated = amble_channel_estimation_interpolated_magnitude.*exp(1i*amble_channel_estimation_interpolated_angle);            

            % create the channel equalization matrix
            channel_estimation = repmat(amble_channel_estimation_interpolated(:,1),1,obj.n_symbols*2);          
            
        end       
        
    elseif strcmp(obj.amble_pos,'pre_end')
        
        % generate a concatenated version of all amble symbols
        Y = [obj.d_k_m_amble_indices.first.pos_f; ...
                obj.d_k_m_amble_indices.second.pos_f; ...
                obj.d_k_m_amble_indices.third.pos_f; ...
                obj.d_k_m_amble_indices.fourth.pos_f];
        X = [obj.d_k_m_amble_indices.first.pos_t; ...
                obj.d_k_m_amble_indices.second.pos_t; ...
                obj.d_k_m_amble_indices.third.pos_t; ...
                obj.d_k_m_amble_indices.fourth.pos_t];
        V = [amble_ch.amb_1_zf; ...
                amble_ch.amb_2_zf; ...
                amble_ch.amb_3_zf; ...
                amble_ch.amb_4_zf];
        
        % interpolate real and imag part
        if strcmp(obj.interpol_ri_ma, 'ri')

            % create interpolant
            F = scatteredInterpolant(X,Y,V, 'linear','nearest');

            % create an interpolation at each point of the frame
            channel_estimation = F(obj.d_k_m_indices.meshgrid.mat_t, obj.d_k_m_indices.meshgrid.mat_f);        

        elseif strcmp(obj.interpol_ri_ma, 'ma')

            % create interpolant
            F_magnitude = scatteredInterpolant(X,Y,abs(V), 'linear','nearest');
            F_angle = scatteredInterpolant(X,Y,V, 'linear','nearest');

            % create an interpolation at each point of the frame
            channel_estimation_magnitude = F_magnitude(obj.d_k_m_indices.meshgrid.mat_t, obj.d_k_m_indices.meshgrid.mat_f);
            channel_estimation_angle = F_angle(obj.d_k_m_indices.meshgrid.mat_t, obj.d_k_m_indices.meshgrid.mat_f);
            channel_estimation_angle = angle(channel_estimation_angle);

            % combine both
            channel_estimation = channel_estimation_magnitude.*exp(1i*channel_estimation_angle);

        else
            error('Unknown interpolation technique.');
        end
    end
    
    % DEBUGGING
    if 1==0
        [xx,yy] = meshgrid(1:1:obj.n_symbols*2,1:1:obj.n_subc);
        figure()
        plot3(xx,yy,real(channel_estimation),'b')
        hold on
        plot3(xx,yy,imag(channel_estimation),'r')
        figure()
        plot3(xx,yy,angle(channel_estimation),'k')
        xlabel('Time Index');
        ylabel('Subcarrier Index');
    end
end