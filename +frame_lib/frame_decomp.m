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

% frame -> d_k_m -> transmitted_bits -> coded_bits
function [coded_bits_rx, transmitted_bits_rx, d_k_m_rx] = frame_decomp(obj)

    % extract d_k_m_rx
    d_k_m_rx = zeros(obj.n_subc, 2*obj.n_symbols);
    d_k_m_rx(:,1:2:end) = real(obj.frame_rx);
    d_k_m_rx(:,2:2:end) = imag(obj.frame_rx);
    
    % extract symbols
	data_symbols_zig_zag = d_k_m_rx(obj.d_k_m_data_indices.lin);
    
    % we can use complex of real modulation, it should be the same in both cases
    if strcmp(obj.modulation_type, 'QAM')
        
        % extract relevant symbols
        if obj.M > 2
            data_symbols_real = data_symbols_zig_zag(1:2:end);
            data_symbols_imag = data_symbols_zig_zag(2:2:end);
        else
            data_symbols_real = data_symbols_zig_zag;
            data_symbols_imag = 0;
        end
        
        % combine real and imaginary part
        data_symbols = data_symbols_real + 1i*data_symbols_imag;        

        % demodulate symbols
        transmitted_bits_rx = qamdemod(data_symbols, obj.M, 'gray', 'OutputType', 'bit','UnitAveragePower', true);        
        
    elseif strcmp(obj.modulation_type, 'PAM')
        
        % for matlabs pamdemod we need to rescale our received symbols
        data_symbols_zig_zag = data_symbols_zig_zag/obj.constell_pam_norm_fac;
        
        % demodulate with matlab
        if obj.M > 2
            transmitted_bits_rx_as_integer = pamdemod(data_symbols_zig_zag, sqrt(obj.M), 0, 'gray');
        else
            transmitted_bits_rx_as_integer = pamdemod(data_symbols_zig_zag, obj.M, 0, 'gray');
        end
        
        % convert integers to bits
        transmitted_bits_rx = de2bi(transmitted_bits_rx_as_integer)';
        
        % reshape bits to one long row vector
        transmitted_bits_rx = reshape(transmitted_bits_rx,[],1);    
        
    else
        
        % should never happen
        strcmp(obj.modulation_type, 'Unknown modulation type.')        
    end
    
    % SECURITY
    if numel(transmitted_bits_rx) ~= obj.n_bits_tx
        error('Wrong number of transmitted bits received.');
    end
    
    % remove additional zeros -> receiver knows from the MAC-header which coding and which block size was used
    coded_bits_rx = transmitted_bits_rx(1:end-obj.n_bits_coded_zeros);
end

