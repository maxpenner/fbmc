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

function [coded_bits] = encode(obj)

    % generate an array containing the blocks
    blocks = zeros(obj.n_bits_user_blk + obj.n_bits_user_zeros_blk, obj.n_blocks);
    
    % write user bits into block array, the tailbits are added here
    blocks(1:obj.n_bits_user_blk, 1:obj.n_blocks) = reshape(obj.user_bits, obj.n_bits_user_blk, obj.n_blocks);     

    % before actual encoding we might have to add some zeros
    coded_bits = zeros(obj.n_blocks*obj.n_bits_coded_blk,1);

    % internal block structure depends on coding
    if strcmp(obj.encoding,'convolutional_7') || strcmp(obj.encoding,'convolutional_9')
        
        % go over each block and encode individually
        for i=1:1:obj.n_blocks
            
            % write start and stop
            write_start = 1 + (i-1)*obj.n_bits_coded_blk;
            write_end = i*obj.n_bits_coded_blk;            
            
            % encode single block of codable bits
            coded_bits_not_interleaved = obj.convEncoder(blocks(:,i));
            
            % interleave
            coded_bits_interleaved = zeros(size(coded_bits_not_interleaved));
            coded_bits_interleaved(obj.block_interl_pattern) = coded_bits_not_interleaved;
            
            % interleave and write to output
            coded_bits(write_start:write_end) = coded_bits_interleaved;
        end
    else
        error('Unknown encoding.');
    end

    % SECURITY
    if numel(coded_bits) ~= obj.n_blocks*obj.n_bits_coded_blk
        error('Incorrect number of coded bits.');
    end
end