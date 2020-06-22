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

% formula (30) from Thein Eurasip

function sto_ = estim_sto(fsync_param, inp_vecs)

    % sto estimation needs amble
    if strcmp(fsync_param.amble_pos, 'pre')
        amble_vec = fsync_param.amble(fsync_param.pilot_idx);
    %elseif strcmp(obj.amble_pos, 'end')
        % TODO
    %elseif strcmp(obj.amble_pos, 'pre_end')
        % TODO
    else
        error('fsync_lib->estim_sto: Amble only defined for PREAMBLE');
    end

    % split for readability
    vec_i_0 = inp_vecs(:,1);
    vec_i_2 = inp_vecs(:,2);
    
    % sto unit is samples
    delta_k = 2;  
    vec_i_0_conj_delta_k = conj(vec_i_0);
    vec_i_2_conj_delta_k = conj(vec_i_2);   
        
    % numerator and denumerator form (30)
    % multiplication is not possible for the last pilot
    numerator_i_0 = vec_i_0(1:end-1).*vec_i_0_conj_delta_k(2:end);
    numerator_i_2 = vec_i_2(1:end-1).*vec_i_2_conj_delta_k(2:end);
    denumerator_i_0_2 = amble_vec(1:end-1).*amble_vec(2:end);
    
    % final sto
    sto_ = sum(numerator_i_0./denumerator_i_0_2) + sum(numerator_i_2./denumerator_i_0_2);
    sto_ = angle(sto_);
    sto_ = fsync_param.n_subc/(2*pi*delta_k)*sto_;
    
    % only integer values are allowed
    if sto_ >= 0
        sto_ = floor(sto_);
    else
        sto_ = ceil(sto_);
    end
end