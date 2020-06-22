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

function [signal_power, error_power, sinr_dB] = points_complex_sinr_db(obj, complex_point_type)

    if strcmp(obj.modulation_type, 'QAM') == false
        error('Complex points SINR can be determined only is QAM is used.');
    end
    
    if strcmp(complex_point_type,'data') == false
        error('Complex points SINR can be determined only for data points.');
    end

    % complex data points we send
    data_symbols_tx = obj.frame(obj.frame_data_indices.lin);

    % what is the average amplitude per complex point -> we must have an average power of 1.0
    amplitude_per_complex_points = 1.0; 

    % from the average amplitude we can calculate the expected average power
    signal_power = amplitude_per_complex_points^2;

    % SECURITY
    % check if power of send complex point is correct
    power_per_complex_point_measured = rms(abs(data_symbols_tx))^2;
    power_per_complex_point_error = (1 - signal_power/power_per_complex_point_measured);
    if power_per_complex_point_error > 5/100
        warning('ON');
        warning('Large offset between set power and measured power for complex point: %.1f %%', power_per_complex_point_error);
        warning('OFF');
    end        

    % complex data points we received
    data_symbols_rx = obj.frame_rx(obj.frame_data_indices.lin);

    % determine the error distance between transmitted complex points and received complex points
    error_euclid_distance = abs(data_symbols_rx - data_symbols_tx);

    % get error power = rms^2 = Effektivwert^2
    error_power = rms(reshape(error_euclid_distance,1,[]))^2;

    % sinr in dB
    sinr_dB = pow2db(signal_power/error_power);       
end

