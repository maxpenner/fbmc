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

function [user_bits_rx] = decode(obj)

    % reshape to coded bits per block
    blocks_coded = reshape(obj.coded_bits_rx, obj.n_bits_coded_blk, obj.n_blocks);
    
    % create an array for the user bits
    user_bits_rx = zeros(obj.n_bits_user_blk*obj.n_blocks, 1);

    % internal block structure depends on coding
    if strcmp(obj.encoding,'convolutional_7') || strcmp(obj.encoding,'convolutional_9')
        
        % decode each block individually
        for i=1:1:obj.n_blocks
            
            % deinterleave block of coded bits
            blocks_coded_deinterleaved = blocks_coded(obj.block_interl_pattern,i);
            
            % decode block
            codable_bits_blk_tmp = obj.vitDecoder(blocks_coded_deinterleaved);
            
            % we know that we only transmit full bytes, so remove possible zeros
            codable_bits_blk_tmp = codable_bits_blk_tmp(1:end-obj.n_bits_user_zeros_blk);
            
            % write to final array
            write_start = 1 + (i-1)*obj.n_bits_user_blk;
            write_end = i*obj.n_bits_user_blk;
            user_bits_rx(write_start:write_end) = codable_bits_blk_tmp;
        end
    else
        error('Unknown decoding.');
    end
    
	% SECURITY
    if numel(user_bits_rx) ~= obj.n_blocks*obj.n_bits_user_blk
        error('Incorrect number of user bits.');
    end
end