% OQAM-OFDM/FBMC for FeelMaTyc
% Author: Maxim Penner (maxim.penner@ikt.uni-hannover.de)
% Date: 2017-05-01
% svn-repo: https://192.168.1.4/svn/feelmatyc/
% copyright (c) 2017 Maxim Penner

% coded_bits -> transmitted_bits -> d_k_m -> frame
function [frame, d_k_m, transmitted_bits] = frame_comp(obj)

    % first append random bits to fill the frame
    transmitted_bits = [obj.coded_bits; randi([0 1],obj.n_bits_coded_zeros,1)];
    
    % SECURITY
    if numel(transmitted_bits) ~= obj.n_bits_tx
        error('Incorrect number of transmitted bits.');
    end
    
    % generate a frame seperated into real and imag part
    d_k_m = zeros(obj.n_subc, 2*obj.n_symbols);
    
    % we can use complex or real modulation
    if strcmp(obj.modulation_type, 'QAM')
        
        % complex symbol generation from transmitted bits
        data_symbols = qammod(transmitted_bits, obj.M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);

        % write the data symbols into d_k_m
        if obj.M > 2
            data_symbols_zig_zag            = zeros(2*numel(data_symbols),1);
            data_symbols_zig_zag(1:2:end)   = real(data_symbols);
            data_symbols_zig_zag(2:2:end)   = imag(data_symbols);
        else
            data_symbols_zig_zag            = real(data_symbols);      
        end        
        
    elseif strcmp(obj.modulation_type, 'PAM')
        
        % group bits into index numbers for constellation diagramm
        if obj.M > 2
            transmitted_bitsgrouped = reshape(transmitted_bits,  log2(obj.M)/2, []);
            bits2num_array = bi2de(reshape(transmitted_bitsgrouped.', [], log2(obj.M)/2));
        else
            transmitted_bitsgrouped = reshape(transmitted_bits,  log2(obj.M), []);
            bits2num_array = bi2de(reshape(transmitted_bitsgrouped.', [], log2(obj.M)));
        end

        % extract the correct constellation point for each integer
        data_symbols_zig_zag = obj.constell_pam_norm(bits2num_array + 1);
        
    else
        strcmp(obj.modulation_type, 'Unknown modulation type.');
    end
    
    % write data symbols
    d_k_m(obj.d_k_m_data_indices.lin) = data_symbols_zig_zag;
    
    % apply warping if wanted
    if strcmp(obj.warp_tx,'none') == true
        % do nothing, just skip
    elseif strcmp(obj.warp_tx,'warp') == true        
        d_k_m = d_k_m.*obj.warp_fac;
    end    
    
    % write pilot and aux
    if numel(obj.pilots) > 0        
        
        % pilot
        d_k_m(obj.d_k_m_pilot_indices.lin) = obj.pilots_val;        
        
        % aux
        d_k_m = equalizer_lib.insertion_aux(obj, d_k_m, obj.pilots(1).f_aux_offset, ...
                                                        obj.pilots(1).t_aux_offset, ...
                                                        obj.d_k_m_pilot_indices.pos_f, ...
                                                        obj.d_k_m_pilot_indices.pos_t);
    end
    
    % write amble and aux
    if obj.amble_len > 0
        
        if strcmp(obj.amble_pos,'pre')
            
            % amble
            d_k_m(obj.d_k_m_amble_indices.first.lin)	= obj.amble_std(obj.d_k_m_amble_indices.first.lin);
            d_k_m(obj.d_k_m_amble_indices.second.lin)	= obj.amble_std(obj.d_k_m_amble_indices.first.lin);
            
            % aux
            d_k_m = equalizer_lib.insertion_aux(obj, d_k_m, 0, 1, obj.d_k_m_amble_indices.second.pos_f, obj.d_k_m_amble_indices.second.pos_t);
            
        elseif strcmp(obj.amble_pos,'end')
            
            % amble
            d_k_m(obj.d_k_m_amble_indices.third.lin)	= obj.amble_std(obj.d_k_m_amble_indices.third.lin_relative);
            d_k_m(obj.d_k_m_amble_indices.fourth.lin)	= obj.amble_std(obj.d_k_m_amble_indices.third.lin_relative);
            
            % aux
            d_k_m = equalizer_lib.insertion_aux(obj, d_k_m, 0, -1, obj.d_k_m_amble_indices.third.pos_f, obj.d_k_m_amble_indices.third.pos_t);
            
        elseif strcmp(obj.amble_pos,'pre_end')
            
            % amble
            d_k_m(obj.d_k_m_amble_indices.first.lin)	= obj.amble_std(obj.d_k_m_amble_indices.first.lin);
            d_k_m(obj.d_k_m_amble_indices.second.lin)	= obj.amble_std(obj.d_k_m_amble_indices.first.lin);
            d_k_m(obj.d_k_m_amble_indices.third.lin)	= obj.amble_std(obj.d_k_m_amble_indices.third.lin_relative);
            d_k_m(obj.d_k_m_amble_indices.fourth.lin)	= obj.amble_std(obj.d_k_m_amble_indices.third.lin_relative);
            
            % aux
            d_k_m = equalizer_lib.insertion_aux(obj, d_k_m, 0, 1, obj.d_k_m_amble_indices.second.pos_f, obj.d_k_m_amble_indices.second.pos_t);
            d_k_m = equalizer_lib.insertion_aux(obj, d_k_m, 0, -1, obj.d_k_m_amble_indices.third.pos_f, obj.d_k_m_amble_indices.third.pos_t);
            
        end        
    end
    
%     % DEBUG
%     %d_k_m(:,:) = 1;
%     
%     d_k_m(:,:) = 0;
%     
%     d_k_m(2:10:end,2:10:end) = -1;
%     d_k_m(2:10:end,3:10:end) = 1;
%     d_k_m(2:10:end,4:10:end) = 1;
%     
%     d_k_m(3:10:end,2:10:end) = 1;
%     %d_k_m(3:10:end,3:10:end) = 1;
%     d_k_m(3:10:end,4:10:end) = -1;    
%    
%     d_k_m(4:10:end,2:10:end) = 1;
%     d_k_m(4:10:end,3:10:end) = 1;
%     d_k_m(4:10:end,4:10:end) = 1;
    
    
    % combine d_k_m to complex frame
    frame = d_k_m(:,1:2:end) + 1i*d_k_m(:,2:2:end);
end











