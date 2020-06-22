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

% m and k are the offset from each point
function [H_m_k_at_each_point] = H_m_k(obj, channel_handle, m_m_dot, k_k_dot, include_factor)

    % determine ideal channel knowledge
    if strcmp(channel_handle.type, 'awgn')

        if m_m_dot == 0 && k_k_dot == 0
            H_m_k_at_each_point = ones(size(obj.oqam_frame));
        else
           	error('Cannot calculate values off the center for awgn.');
        end

    elseif strcmp(channel_handle.type, 'deterministic')

        H_m_k_at_each_point = ones(size(obj.oqam_frame));
        
        warning('ON');
        warning('Calculating channel coefficient H_k_m not implemented for deterministic channel.');
        warning('ON');

    elseif strcmp(channel_handle.type, 'statistical')     
        
        % slider
            
        % define a perfect slider
        perfect_slider_statistical = [];
        
        % we need the factor K
        K = obj.n_subc;
        T = obj.n_subc;

        % we define filter as p in equations
        p = obj.filter_func;
        n_p = numel(p);

        % how many path delays do we have?
        n_pathDelays = size(channel_handle.s_rayleigh_gains, 1);

        % init the slider p*p_delayed
        perfect_slider_statistical.p_p_delayed = repmat(p, n_pathDelays, 1);
        for i = 1:1:n_pathDelays
            
            % delay due to tau
            tau = i-1;
            
            % delay for specific symbol
            full_delay = tau + m_m_dot*T/2;
            
            % multiply windows
            if full_delay <= -n_p
                
                p_delayed = zeros(1, n_p);
                
            elseif full_delay > -n_p && full_delay < 0
                
                p_delayed = [p(-full_delay+1:end), zeros(1, -full_delay)];
                
            elseif full_delay >= 0 && full_delay <= n_p
                
                p_delayed = [zeros(1, full_delay), p(1:end-full_delay)];
                
            elseif full_delay > n_p
                
                p_delayed = zeros(1, n_p);
                
            end
            
            % multiply windows
            perfect_slider_statistical.p_p_delayed(i,:) = perfect_slider_statistical.p_p_delayed(i,:).*p_delayed;
        end

        % init the slider to the exponential function with k_dot * tau
        k_dot_vec = -K/2 : 1 : (K/2-1);
        kappa_plus_k_dot_vec = k_dot_vec + k_k_dot;
        kappa_plus_k_dot_vec = repmat(kappa_plus_k_dot_vec', 1, n_pathDelays);

        % create a vector of delays
        tau_mat = 0:1:(n_pathDelays-1);
        tau_mat = repmat(tau_mat, K, 1);

        % create the exponent of the e-function
        if strcmp(include_factor,'include') == true
            
            exponent = 1i*2*pi/K*kappa_plus_k_dot_vec.*(-tau_mat - m_m_dot*T/2);
            
        elseif strcmp(include_factor,'include not') == true
            
            exponent = 1i*2*pi/K*kappa_plus_k_dot_vec.*(-tau_mat - 0*m_m_dot*T/2);
            
        else
            error('Unknown factor behaviour.');
        end

        % create the second slider
        perfect_slider_statistical.e_slider = exp(exponent);
        
        %% use the slider and slide over all symbols

        % initialize the inv channel to zeros
        H_m_k_at_each_point = zeros(size(obj.oqam_frame));

        % determine how many symbols we have
        n_m_dot = size(obj.oqam_frame, 2);

        % next we have to calculate the multiplication of p*p_delayed and the channel coefficients
        for m_dot = 0:1:(n_m_dot-1)

            % what's the matlab index
            m_dot_matlab = m_dot + 1;
            
            % we integrate over n strich
            n_strich = 0 : 1 : (numel(obj.filter_func)-1);
            
            % mix up the multiplication of the filter
            p_p_delayed_mixed = perfect_slider_statistical.p_p_delayed.*exp(1i*2*pi/K*k_k_dot*n_strich);
            
            % h*p*p_delayed
            gains_range_start = 1+m_dot*T/2;
            gains_range_end = gains_range_start + numel(obj.filter_func)-1;            
            
            % multiply with channel coefficient
            summation_over_n = channel_handle.s_rayleigh_gains(:, gains_range_start:gains_range_end).*p_p_delayed_mixed;

            % sum over n
            summation_over_n = sum(summation_over_n,2);
            summation_over_n = summation_over_n.';
            summation_over_n = repmat(summation_over_n, K, 1);

            % mutliply with the e-function slider
            H_k_dot_m_dot_00_over_each_subcarrier = perfect_slider_statistical.e_slider.*summation_over_n;
            H_k_dot_m_dot_00_over_each_subcarrier = sum(H_k_dot_m_dot_00_over_each_subcarrier,2);

            % write into inv_channel_estimation
            H_m_k_at_each_point(:,m_dot_matlab) = H_k_dot_m_dot_00_over_each_subcarrier;
        end
    else
        error('Unknown channel type.');
    end
    
	% The channel object reduces the amplitude of the incoming samples.
    % After this reduction, the channel preserves the power.
    % So the amplitude change is not included in the channel and has to be added here.
    if channel_handle.amp ~= 1.0
        H_m_k_at_each_point = H_m_k_at_each_point*channel_handle.amp;
    end
end

