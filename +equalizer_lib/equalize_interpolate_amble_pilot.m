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

function [channel_estimation] = equalize_interpolate_amble_pilot(obj, amble_ch, pilot_ch)

    % channel estimation depends on our frame structure
    if strcmp(obj.amble_pos,'pre')
        
        % generate a concatenated version of all amble symbols
        Y = [obj.d_k_m_amble_indices.first.pos_f; ...
                obj.d_k_m_amble_indices.second.pos_f;];
        X = [obj.d_k_m_amble_indices.first.pos_t; ...
                obj.d_k_m_amble_indices.second.pos_t];
        V = [amble_ch.amb_1_zf; ...
                amble_ch.amb_2_zf];

    elseif strcmp(obj.amble_pos,'end')

        % generate a concatenated version of all amble symbols
        Y = [obj.d_k_m_amble_indices.third.pos_f; ...
                obj.d_k_m_amble_indices.fourth.pos_f;];
        X = [obj.d_k_m_amble_indices.third.pos_t; ...
                obj.d_k_m_amble_indices.fourth.pos_t];
        V = [amble_ch.amb_3_zf; ...
                amble_ch.amb_4_zf];        
        
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

    end

    % append the channel estimation points from our pilot symbols as well
    Y = [Y; obj.d_k_m_pilot_indices.pos_f];
    X = [X; obj.d_k_m_pilot_indices.pos_t];
    V = [V; pilot_ch];
    
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