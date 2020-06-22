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

function [n_bits_tx, n_bits_coded_zeros, n_blocks, n_bits_coded_blk, n_bits_codable_blk_tmp, n_bits_user_zeros_blk, n_bits_user_blk] = bits_count(obj)
%function [n_bits_tx, n_bits_coded, n_bits_coded_zeros, n_bits_codable, n_bits_user, n_bits_user_zeros, n_blocks] = bits_count(obj)

    % determine the number of bits we can actually transmit
    if obj.M == 2
        n_bits_tx = numel(obj.d_k_m_data_indices.lin)*log2(obj.M);
    else
        n_bits_tx = numel(obj.d_k_m_data_indices.lin)*log2(obj.M)/2;
    end
    
    % first we need to determine a temporary size for our blocks
    if obj.n_bits_codable_blk_wish == inf
        
        % calculate a block size so that after coding we can fit exactly one block in the number of transmitable bits
        n_bits_codable_blk_tmp = floor(n_bits_tx/obj.denumerator)*obj.numerator;
        
    elseif obj.n_bits_codable_blk_wish > 0 && obj.n_bits_codable_blk_wish < inf
        
        % we keep the number we actually programmed
        n_bits_codable_blk_tmp = obj.n_bits_codable_blk_wish;
        
    else
        error('Impossible number of bits per block.');        
    end
    
    % determine the size of a block after coding
    n_bits_coded_blk = floor(n_bits_codable_blk_tmp*obj.denumerator/obj.numerator);
    
    % SECURITY
    if n_bits_coded_blk ~= n_bits_codable_blk_tmp*obj.denumerator/obj.numerator
        error('The number of codable bits must be changed for this code rate.');
    end    
    
    % determine how many frames fit in the number of transmit bits
    n_blocks = floor(n_bits_tx/n_bits_coded_blk);
    
    % SECURITY
    if n_blocks <= 0
        error('Number of blocks is too low.');
    end
    
    % SECURITY
    if obj.n_bits_codable_blk_wish == inf && n_blocks ~= 1
        error('With infinite block size exactly one block must be generated.');
    end
    
    % now we can calculate how many bits have to be appended after coding to ALL blocks
    n_bits_coded_zeros = n_bits_tx - n_blocks*n_bits_coded_blk;
      
    % next we have to determine how many user bits fit in one block
    if strcmp(obj.encoding,'convolutional_7')
        
        % 6 bits of the user bits need to be zeros for the viterbi decoder.
        % Also, we only transmit bytes.
        n_bits_user_blk = floor((n_bits_codable_blk_tmp - 6)/8)*8;
        
        % we add at least 6 zeros to the user bits, but it might be even more
        n_bits_user_zeros_blk = n_bits_codable_blk_tmp - n_bits_user_blk;
        
        % SECURITY
        if n_bits_user_zeros_blk < 6
            error('Not enough zeros appended.');
        end
        
    elseif strcmp(obj.encoding,'convolutional_9')
        
        % 8 bits of the user bits need to be zeros for the viterbi decoder.
        % Also, we only transmit bytes.
        n_bits_user_blk = floor((n_bits_codable_blk_tmp - 8)/8)*8;
        
        % we add at least 6 zeros to the user bits, but it might be even more
        n_bits_user_zeros_blk = n_bits_codable_blk_tmp - n_bits_user_blk;
        
        % SECURITY
        if n_bits_user_zeros_blk < 8
            error('Not enough zeros appended.');
        end
        
    else                
        error('Undefined encoding.');
    end
    
    % SECURITY
    if obj.n_bits_codable_blk_wish > 0 && obj.n_bits_codable_blk_wish < inf
        if obj.n_bits_codable_blk_wish ~= n_bits_codable_blk_tmp
            error('Number of codable bits has changed.');
        end
    end    
    
    % SECURITY
    if n_bits_codable_blk_tmp ~= n_bits_user_blk + n_bits_user_zeros_blk
        error('Codable bit number is incorrect.');
    end
    if n_bits_coded_blk ~= n_bits_codable_blk_tmp*obj.denumerator/obj.numerator
        error('Coded bit number is incorrect.');
    end
    if n_bits_tx ~= n_bits_coded_blk*n_blocks + n_bits_coded_zeros
        error('Transmit bit number is incorrect.');
    end
    
%     % next we have to determine how many user bits fit in one block
%     if strcmp(obj.encoding,'none')
%         
%         % we only transmit full bytes
%         n_bits_user = floor(n_bits_tx/8)*8;
%         
%         % the other numbers are
%         n_bits_coded = n_bits_user;        
%         n_bits_coded_zeros = n_bits_tx - n_bits_coded;
%         n_bits_codable = n_bits_user;
%         n_bits_user_zeros = n_bits_codable - n_bits_user;
%         
%     elseif strcmp(obj.encoding,'convolutional_7')
%         
%         % we need a multiple of the encoder output
%         n_bits_coded = floor(n_bits_tx/obj.denumerator)*obj.denumerator;
%         
%         % we might need to append a few zeros
%         n_bits_coded_zeros = n_bits_tx - n_bits_coded;
%         
%         % determine number of codable bits
%         n_bits_codable = n_bits_coded/obj.denumerator*obj.numerator;
%         
%         % 6 bits of the user bits need to be zeros for the viterbi decoder.
%         % Also, we only transmit bytes.
%         n_bits_user = floor((n_bits_codable - 6)/8)*8;
%         
%         % we add at least 6 zeros to the user bits, but it might be even more
%         n_bits_user_zeros = n_bits_codable - n_bits_user;
%         
%         % SECURITY
%         if n_bits_user_zeros < 6
%             error('Not enough zeros appended.');
%         end
%         
%     elseif strcmp(obj.encoding,'blockcode')
%         
%         % TODO
%         
%     else
%         error('Undefined encoding.');
%     end
%     
%     % SECURITY
%     if n_bits_tx ~= n_bits_coded + n_bits_coded_zeros
%         error('Transmit bit number is incorrect.');
%     end
%     if n_bits_coded ~= n_bits_codable*obj.denumerator/obj.numerator
%         error('Coded bit number is incorrect.');
%     end
%     if n_bits_codable ~= n_bits_user + n_bits_user_zeros
%         error('Codable bit number is incorrect.');
%     end
end
