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

function [pilots_val, pilots_val_pattern] = pilot_values(obj)

    % generate a random stream of 1 and -1 which always starts equally
    pilots_stream = randi(RandStream('mt19937ar','Seed',12345), [0 1], numel(obj.d_k_m_pilot_indices.lin), 1)*2-1;

    % each pilot is real in d_k_m frame with an average power of 0.5
    pilot_amplitude_unboosted = 1/sqrt(2);
    
    % convert stream of 1 and -1 to correct pilot amplitudes
    pilot_unboosted = pilot_amplitude_unboosted*pilots_stream;

    % SECURITY
    % check if pilots have an average power of 0.5
    power_pilots_unboosted = (rms(pilot_unboosted))^2;
    if abs((1 - power_pilots_unboosted/0.5)) > 1e-3
        error('Unboosted pilots must have an average power of 0.5.');
    end                

    % boost pilot amplitude
    pilots_val = pilot_unboosted*obj.pilots_amplitude_boost;

    % calculate a pilot version with real-imag-pattern
    pilots_val_pattern = pilots_val.*((1i).^(mod(obj.d_k_m_pilot_indices.pos_f - 1 + obj.d_k_m_pilot_indices.pos_t - 1, 2)));

    % SECURITY
    % check if pilots have correct boosted amplitude
    pilot_avg_amplitude = mean(abs(pilots_val_pattern));
    if abs((1 - pilot_avg_amplitude/(pilot_amplitude_unboosted*obj.pilots_amplitude_boost))) > 1e-3
        error('Boosted pilots with pattern have incorrect average amplitude.');
    end
end