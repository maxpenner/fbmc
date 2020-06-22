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

function [frame_data_indices, d_k_m_indices, d_k_m_data_indices, d_k_m_amble_indices, d_k_m_pilot_indices, d_k_m_statistics] = patterns(obj)
    
    %% we need to determine the position of these complex points
    shape_data_zero     = 1;
    shape_data          = 2;
    shape_amble_zero    = 3;
    shape_amble         = 4;
    shape_amble_aux   	= 5;
    shape_guard         = 6;
    shape_hole          = 7;
    shape_pilot_zero    = 8;
    shape_pilot         = 9;
    shape_pilot_aux     = 10;

    %% generate a shape of the entire d_k_m = complex frame split into real and imag
    d_k_m_shape = zeros(obj.n_subc, 2*obj.n_symbols);
    d_k_m_shape_dim_freq = obj.n_subc;
    d_k_m_shape_dim_time = 2*obj.n_symbols;
    d_k_m_shape_dim = [d_k_m_shape_dim_freq, d_k_m_shape_dim_time];
    
    %% assume all d_k_m points are data zeros
    d_k_m_shape(:,:) = shape_data_zero;

    %% assign complex data points
    % with BPSK we are only transmitting in the real part, so we are taking only each second vector
    if obj.M == 2
        d_k_m_shape(:,1:2:end) = shape_data;
    else
        d_k_m_shape(:,:) = shape_data;
    end
    
    %% overwrite with points, where we deactivate bits due to warping
    if strcmp(obj.warp_tx,'warp') == true
        
        % SECURITY
        if numel(obj.warp_fac) == 0
            error('Trying to warp with no data.');
        end
        
        % find all indices, where warp_fac is zero
        [warp_zero_pos_f, warp_zero_pos_t, ~] = find(obj.warp_fac == 0);
        
        % convert to indices
        warp_zero_indices = sub2ind([d_k_m_shape_dim_freq, d_k_m_shape_dim_time], warp_zero_pos_f, warp_zero_pos_t);
        
        % write data zero where warping is zero
        d_k_m_shape(warp_zero_indices) = shape_data_zero;
    end

    %% next write amble complex points and amble aux
    amble_offset_left = 0;
    amble_offset_right = 0;
    if obj.amble_len > 0
        if strcmp(obj.amble_pos, 'pre')
            d_k_m_shape(:,1:2*obj.amble_len)            = shape_amble_zero;
            d_k_m_shape(1:2:end,1)                      = shape_amble;
            d_k_m_shape(1:2:end,3)                      = shape_amble;
            d_k_m_shape(1:2:end,4)                      = shape_amble_aux;
            amble_offset_left                           = 2*obj.amble_len;
            amble_offset_right                          = 0;
        elseif strcmp(obj.amble_pos, 'end')
            d_k_m_shape(:,end-2*obj.amble_len+1:end)	= shape_amble_zero;
            d_k_m_shape(1:2:end,end)                    = shape_amble;
            d_k_m_shape(1:2:end,end-2)                  = shape_amble;
            d_k_m_shape(1:2:end,end-3)                  = shape_amble_aux;
            amble_offset_left                           = 0;
            amble_offset_right                          = 2*obj.amble_len;            
        elseif strcmp(obj.amble_pos, 'pre_end')
            d_k_m_shape(:,1:2*obj.amble_len)            = shape_amble_zero;
            d_k_m_shape(1:2:end,1)                      = shape_amble;
            d_k_m_shape(1:2:end,3)                      = shape_amble;
            d_k_m_shape(1:2:end,4)                      = shape_amble_aux;
            d_k_m_shape(:,end-2*obj.amble_len+1:end)	= shape_amble_zero;
            d_k_m_shape(1:2:end,end)                    = shape_amble;
            d_k_m_shape(1:2:end,end-2)                  = shape_amble;
            d_k_m_shape(1:2:end,end-3)                  = shape_amble_aux;
            amble_offset_left                           = 2*obj.amble_len;
            amble_offset_right                          = 2*obj.amble_len;
        else
            error('Unknown amble position.');
        end
    end
    
    %% next write the guard positions and the DC
    % we can have a frame without guards  
    d_k_m_shape(1:obj.n_guard_l, :) = shape_guard;          % left spectrum edge
    d_k_m_shape(end-obj.n_guard_r+1:end, :)	= shape_guard;  % right spectrum edge
    
	% we can also deactivate the DC carrier
    if obj.deac_DC == true
        d_k_m_shape(obj.n_subc/2+1, :) = shape_guard;
    end
    
    %% write all holes for all schemes
    for i=1:1:numel(obj.holes)
        hole_time_start	= 1 + 2*(obj.holes(i).time_start-1);    % transform frame to d_k_m index (1 frame -> 3 d_k_m)
        hole_time_end   = 2 + 2*(obj.holes(i).time_end-1);      % transform frame to d_k_m index
        hole_freq_start = obj.holes(i).freq_start;
        hole_freq_end   = obj.holes(i).freq_end;
        
        % put the whole in the spectrum
        d_k_m_shape(hole_freq_start:hole_freq_end, hole_time_start:hole_time_end) = shape_hole;
    end    
    
    %% next write embedded pilot and aux pilot
    % SECURITY
    % all pilot schemes must use the same surrounding aux
    if numel(obj.pilots) > 0
        for i=2:1:numel(obj.pilots)
            if obj.pilots(1).aux_indices ~= obj.pilots(i).aux_indices
                error('All pilot schemes must use the same aux indices.');
            end
        end
    end
    
    % we can define multiple schemes
    for i=1:1:numel(obj.pilots)
        
        % time
        for time_pos = 1 + amble_offset_left + 2*obj.pilots(i).time_offset : 2*obj.pilots(i).time_step : d_k_m_shape_dim_time - amble_offset_right
            
            % freq
            for freq_pos = 1 + obj.n_guard_l + obj.pilots(i).freq_offset : obj.pilots(i).freq_step : d_k_m_shape_dim_freq - obj.n_guard_r
                
                % set pilot shape only at data shapes
                if d_k_m_shape(freq_pos, time_pos) == shape_data || d_k_m_shape(freq_pos, time_pos) == shape_data_zero
                    d_k_m_shape(freq_pos, time_pos) = shape_pilot;
                else
                    continue;
                end
                
                % set auxiliary pilots
                for j = 1:1:numel(obj.pilots(i).t_aux_offset)
                    freq_pos_aux = freq_pos + obj.pilots(i).f_aux_offset(j);
                    time_pos_aux = time_pos + obj.pilots(i).t_aux_offset(j);
                    
                    % try setting aux pilot
                    if d_k_m_shape(freq_pos_aux, time_pos_aux) == shape_data || d_k_m_shape(freq_pos_aux, time_pos_aux) == shape_data_zero
                        d_k_m_shape(freq_pos_aux, time_pos_aux) = shape_pilot_aux;
                    else
                        error('Trying to put a aux shape on a non-data shape.');
                    end
                end
                
                % set blockings
                for j = 1:1:numel(obj.pilots(i).block_t_offset)
                    freq_pos_block = freq_pos + obj.pilots(i).block_f_offset(j);
                    time_pos_block = time_pos + obj.pilots(i).block_t_offset(j);
                    
                    % try setting aux pilot
                    if d_k_m_shape(freq_pos_block, time_pos_block) == shape_data || d_k_m_shape(freq_pos_block, time_pos_block) == shape_data_zero
                        d_k_m_shape(freq_pos_block, time_pos_block) = shape_pilot_zero;
                    else
                        error('Trying to put a pilot blocking shape on a non-data shape.');
                    end
                end
            end
        end
    end
    
    %% create indices for each complex point
    [d_k_m_indices.pos_f, d_k_m_indices.pos_t] = ind2sub(d_k_m_shape_dim, 1:1:d_k_m_shape_dim_freq*d_k_m_shape_dim_time);
    [d_k_m_indices.meshgrid.mat_t, d_k_m_indices.meshgrid.mat_f] = meshgrid(1:1:d_k_m_shape_dim_time, 1:1:d_k_m_shape_dim_freq);
    
    %% determine indices of data shape
    % determine number of data points in d_k_m
    n_data_points_d_k_m = numel(find(d_k_m_shape - shape_data == 0));
    
    % extract the indices within frame if real OR imag part was used
    frame_data_indices.lin = zeros(n_data_points_d_k_m,1);
    
    % preallocate memory for d_k_m
    d_k_m_data_indices.lin = zeros(n_data_points_d_k_m,1);
    
    % now iterate over each real vector (1, 3, 5...) and collect indices real and imag vectors
    idx = 0;
    idx_frame = 0;
    for i=1:2:d_k_m_shape_dim_time
        for j=1:1:d_k_m_shape_dim_freq
            
            % check complex point in real vector
            if d_k_m_shape(j,i) == shape_data
                idx = idx + 1;
                d_k_m_data_indices.lin(idx) = (i-1)*d_k_m_shape_dim_freq + j;               
            end
            
            % check complex point in imag vector right next to it
            if d_k_m_shape(j,i+1) == shape_data
                idx = idx + 1;
                d_k_m_data_indices.lin(idx) = i*d_k_m_shape_dim_freq + j;
            end
            
            % check complex point in real vector
            if d_k_m_shape(j,i) == shape_data || d_k_m_shape(j,i+1) == shape_data
                idx_frame = idx_frame + 1;
                frame_data_indices.lin(idx_frame) = (i-1)/2*d_k_m_shape_dim_freq + j;
            end            
        end
    end
    
    % cut frame indices to relevant part
    frame_data_indices.lin = frame_data_indices.lin(1:idx_frame);
    
    % SECURITY
    if strcmp(obj.modulation_type, 'QAM') == true && numel(frame_data_indices.lin)*2 ~= numel(d_k_m_data_indices.lin)
        error('For QAM frame must contain twice as many indices as d_k_m.');
    end    
    
    % SECURITY
    if idx ~= numel(d_k_m_data_indices.lin)
        error('Collected data points and idx counter are not equal.');
    end    
    
    % SECURITY
    if numel(d_k_m_data_indices.lin) == 0
        error('No complex data points defined.');
    end
    
    % convert linear positions to coordinates
    [frame_data_indices.pos_f, frame_data_indices.pos_t] = ind2sub([d_k_m_shape_dim_freq, d_k_m_shape_dim_time/2], frame_data_indices.lin);
    [d_k_m_data_indices.pos_f, d_k_m_data_indices.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_data_indices.lin);
    
    %% determine indices of amble shape
    d_k_m_amble_indices = [];
    if obj.amble_len > 0
        if strcmp(obj.amble_pos,'pre')
            
          	first = find(d_k_m_shape(:,1) - shape_amble == 0);
            second = find(d_k_m_shape(:,3) - shape_amble == 0);
            
            % SECURITY
            if numel(first) == 0 || numel(second) == 0
                error('No amble indices defined.');
            end
            
            % SECURITY
            if sum(abs(first - second)) ~= 0
                error('First two ambles vectors are not equal.');
            end
            
            % calculate linear indices within d_k_m shape
            d_k_m_amble_indices.first.lin = first;
            d_k_m_amble_indices.second.lin = second + 2*d_k_m_shape_dim_freq;
            
            % convert to absolute position in entire d_k_m shape
            [d_k_m_amble_indices.first.pos_f, d_k_m_amble_indices.first.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.first.lin);
            [d_k_m_amble_indices.second.pos_f, d_k_m_amble_indices.second.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.second.lin);
            
        elseif strcmp(obj.amble_pos,'end')
            
        	third = find(d_k_m_shape(:,end-2) - shape_amble == 0);
            fourth = find(d_k_m_shape(:,end) - shape_amble == 0);
            
            % SECURITY
            if numel(third) == 0 || numel(fourth) == 0
                error('No amble indices defined.');
            end
            
            % SECURITY
            if sum(third - fourth) ~= 0
                error('Last two ambles vectors are not equal.');
            end
            
            % calculate linear indices within d_k_m shape
            d_k_m_amble_indices.third.lin = d_k_m_shape_dim_freq*(d_k_m_shape_dim_time - 3) + third;
            d_k_m_amble_indices.fourth.lin = d_k_m_shape_dim_freq*(d_k_m_shape_dim_time - 1) + fourth;
            
            % convert to absolute position in entire d_k_m shape
            [d_k_m_amble_indices.third.pos_f, d_k_m_amble_indices.third.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.third.lin);
            [d_k_m_amble_indices.fourth.pos_f, d_k_m_amble_indices.fourth.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.fourth.lin);
            
            % add indices relative in the third vector
            d_k_m_amble_indices.third.lin_relative = third;
            
        elseif strcmp(obj.amble_pos,'pre_end')
            
          	first = find(d_k_m_shape(:,1) - shape_amble == 0);
            second = find(d_k_m_shape(:,3) - shape_amble == 0);
        	third = find(d_k_m_shape(:,end-2) - shape_amble == 0);
            fourth = find(d_k_m_shape(:,end) - shape_amble == 0);            
            
            % SECURITY
            if numel(first) == 0 || numel(second) == 0
                error('No amble indices defined.');
            end
            
            % SECURITY
            if sum(abs(first - second)) ~= 0
                error('First two ambles vectors are not equal.');
            end
            
            % SECURITY
            if numel(third) == 0 || numel(fourth) == 0
                error('No amble indices defined.');
            end
            
            % SECURITY
            if sum(third - fourth) ~= 0
                error('Last two ambles vectors are not equal.');
            end            
            
            % calculate linear indices within d_k_m shape
            d_k_m_amble_indices.first.lin = first;
            d_k_m_amble_indices.second.lin = second + 2*d_k_m_shape_dim_freq;
            d_k_m_amble_indices.third.lin = d_k_m_shape_dim_freq*(d_k_m_shape_dim_time - 3) + third;
            d_k_m_amble_indices.fourth.lin = d_k_m_shape_dim_freq*(d_k_m_shape_dim_time - 1) + fourth;            
            
            % convert to absolute position in entire d_k_m shape
            [d_k_m_amble_indices.first.pos_f, d_k_m_amble_indices.first.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.first.lin);
            [d_k_m_amble_indices.second.pos_f, d_k_m_amble_indices.second.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.second.lin);
            [d_k_m_amble_indices.third.pos_f, d_k_m_amble_indices.third.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.third.lin);
            [d_k_m_amble_indices.fourth.pos_f, d_k_m_amble_indices.fourth.pos_t] = ind2sub(d_k_m_shape_dim, d_k_m_amble_indices.fourth.lin);
            
            % add indices relative in the third vector
            d_k_m_amble_indices.third.lin_relative = third;            
        end
    end
    
    %% determine indices of pilot indices    
    d_k_m_pilot_indices = [];
    if numel(obj.pilots) > 0
        
        [d_k_m_pilot_indices.pos_f, d_k_m_pilot_indices.pos_t, ~] = find(d_k_m_shape - shape_pilot == 0);
        
        % SECURITY
        if numel(d_k_m_pilot_indices) == 0
            error('No pilot indices defined.');
        end
        
        % convert coordinates to linear indices
        d_k_m_pilot_indices.lin = sub2ind(d_k_m_shape_dim, d_k_m_pilot_indices.pos_f, d_k_m_pilot_indices.pos_t);
    end
    
    %% calculate some statistics about how d_k_m is used
    d_k_m_statistics.d_k_m_points_total = numel(d_k_m_shape);
    d_k_m_statistics.d_k_m_data_points_total = numel(find(d_k_m_shape == shape_data));
    d_k_m_statistics.d_k_m_amble_points_total = numel(find(d_k_m_shape == shape_amble));
    d_k_m_statistics.d_k_m_amble_aux_points_total = numel(find(d_k_m_shape == shape_amble_aux));
    d_k_m_statistics.d_k_m_pilot_points_total = numel(find(d_k_m_shape == shape_pilot));
    d_k_m_statistics.d_k_m_pilot_aux_points_total = numel(find(d_k_m_shape == shape_pilot_aux));
    
    sum_used_d_k_m_points = d_k_m_statistics.d_k_m_data_points_total + ...
                            d_k_m_statistics.d_k_m_amble_points_total + ...
                            d_k_m_statistics.d_k_m_amble_aux_points_total + ...
                            d_k_m_statistics.d_k_m_pilot_points_total + ...
                            d_k_m_statistics.d_k_m_pilot_aux_points_total;
    
    d_k_m_statistics.usage = sum_used_d_k_m_points/d_k_m_statistics.d_k_m_points_total;    
    
    %% Debugging Plot
    if 1==0
        figure()
        clf()
        imagesc(d_k_m_shape)
        axis image
        axis ij
        colormap jet
    end
end

