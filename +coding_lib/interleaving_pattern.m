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

function [block_interl_pattern] = interleaving_pattern(obj)

    % the interleaving pattern is applied over the coded bits of a block
    if obj.block_interl_active == true
        
        % we interleave over all coded bits of one block
        block_interl_pattern = randperm(obj.n_bits_coded_blk);
        
    else
        
        % we assume a simple pattern of ones
        block_interl_pattern = 1:1:obj.n_bits_coded_blk;
    end
    
    % SECURITY
    if unique(block_interl_pattern) ~= obj.n_bits_codable_blk
        error('Interleaving pattern is not unique.');
    end
end