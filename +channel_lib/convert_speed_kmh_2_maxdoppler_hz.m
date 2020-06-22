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

function [max_doppler_hz] = convert_speed_kmh_2_maxdoppler_hz(speed_kmh,carrier_frequency)

    % convert speed from km/h to m/s
    speed_ms = speed_kmh*1000/(60*60);

    % calculate wavelength from carrier frequency
    wavelength = physconst('LightSpeed')/carrier_frequency;
    
    % calculate the maximum doppler frequency
    max_doppler_hz = speed_ms/wavelength;
end

